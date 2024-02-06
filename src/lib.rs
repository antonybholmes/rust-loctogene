use std::{error::Error, fmt};

use dna::Location;
use r2d2_sqlite::SqliteConnectionManager;

use serde::Serialize;

mod tests;

const WITHIN_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
    FROM genes 
    WHERE level = ? AND chr = ? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) 
    ORDER BY start ASC"#;

const IN_EXON_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
    FROM genes 
    WHERE level=3 AND gene_id=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) 
    ORDER BY start ASC"#;

const IN_PROMOTER_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
    FROM genes 
    WHERE level=2 AND gene_id=? AND chr=? AND ? >= stranded_start - ? AND ? <= stranded_start + ? 
    ORDER BY start ASC"#;

const CLOSEST_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, stranded_start - ? 
	FROM genes
	WHERE level=? AND chr=?
	ORDER BY ABS(stranded_start - ?) 
	LIMIT ?"#;

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Plus = 1,
    Neg = 2,
}

impl From<&str> for Strand {
    fn from(level: &str) -> Self {
        match level {
            "-" => Strand::Plus,
            _ => Strand::Plus,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Neg => write!(f, "-"),
            _ => write!(f, "+"),
        }
    }
}

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Level {
    Gene = 1,
    Transcript = 2,
    Exon = 3,
}

impl From<&str> for Level {
    fn from(level: &str) -> Self {
        match level {
            "transcript" => Level::Transcript,
            "exon" => Level::Exon,
            "2" => Level::Transcript,
            "3" => Level::Exon,
            _ => Level::Gene,
        }
    }
}

impl From<u8> for Level {
    fn from(level: u8) -> Self {
        match level {
            2 => Level::Transcript,
            3 => Level::Exon,
            _ => Level::Gene,
        }
    }
}

impl fmt::Display for Level {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Level::Gene => write!(f, "Gene"),
            Level::Transcript => write!(f, "Transcript"),
            Level::Exon => write!(f, "Exon"),
        }
    }
}

#[derive(Serialize, Debug, PartialEq, Eq, Clone, Copy)]
pub struct TSSRegion {
    pub offset_5p: i32,
    pub offset_3p: i32,
}

impl TSSRegion {
    pub fn new(offset_5p: i32, offset_3p: i32) -> TSSRegion {
        return TSSRegion {
            offset_5p,
            offset_3p,
        };
    }
}

pub const DEFAULT_TSS_REGION: TSSRegion = TSSRegion {
    offset_5p: -2000,
    offset_3p: 1000,
};

#[derive(Serialize, Debug, PartialEq, Eq, Clone)]
pub struct GenomicFeature {
    pub id: u32,
    pub chr: String,
    pub start: i32,
    pub end: i32,
    pub strand: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub dist: i32,
}

// #[derive(Serialize)]
// pub struct GenomicFeatures {
//     pub level: Level,
//     pub features: Vec<GenomicFeature>,
// }

//const NO_FEATURES: [Features; 0] = [] .to_vec();

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

pub struct Loctogene {
    pool: r2d2::Pool<SqliteConnectionManager>,
}

impl Loctogene {
    pub fn new(file: &str) -> Result<Loctogene, Box<dyn Error>> {
        // let db: Connection = match Connection::open(file) {
        //     Ok(db) => db,
        //     Err(err) => return Err(format!("{}", err)),
        // };

        let manager: SqliteConnectionManager = SqliteConnectionManager::file(file);

        let pool: r2d2::Pool<SqliteConnectionManager> = r2d2::Pool::builder().build(manager)?;

        Ok(Loctogene { pool })
    }

    // pub fn get_genes_within_stranded(
    //     &self,
    //     location: &Location,
    //     strand: Strand,
    //     level: Level,
    // ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
    //     let mid: u32 = location.mid();

    //     let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

    //     let mut stmt: rusqlite::Statement<'_> = pool.prepare(WITHIN_GENE_STRANDED_SQL)?;

    //     let mapped_rows = stmt.query_map(
    //         rusqlite::params![
    //             mid,
    //             level as u8,
    //             location.chr,
    //             strand.to_string(),
    //             location.start,
    //             location.start,
    //             location.end,
    //             location.end
    //         ],
    //         |row: &rusqlite::Row<'_>| row_to_feature(row),
    //     )?;

    //     let features: Vec<GenomicFeature> = mapped_rows
    //         .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
    //         .collect::<Vec<GenomicFeature>>();

    //     Ok(features)
    // }

    pub fn get_genes_within(
        &self,
        location: &Location,
        level: Level,
    ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
        let mid: i32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

        let mut stmt: rusqlite::Statement<'_> = pool.prepare(WITHIN_GENE_SQL)?;

        let mapped_rows = stmt.query_map(
            rusqlite::params![
                mid,
                level as u8,
                location.chr,
                location.start,
                location.start,
                location.end,
                location.end
            ],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        )?;

        let features: Vec<GenomicFeature> = mapped_rows
            .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
            .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns the exons that a location is in within a particular gene. Useful
    // for determining if a gene is exonic or not.
    pub fn in_exon(
        &self,
        location: &Location,
        gene_id: &str,
    ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
        let mid: i32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

        let mut stmt: rusqlite::Statement<'_> = pool.prepare(IN_EXON_SQL)?;

        let mapped_rows = stmt.query_map(
            rusqlite::params![
                mid,
                gene_id,
                location.chr,
                location.start,
                location.start,
                location.end,
                location.end
            ],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        )?;

        let features: Vec<GenomicFeature> = mapped_rows
            .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
            .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns a list of features if location is in tss of specific gene
    pub fn in_promoter(
        &self,
        location: &Location,
        gene_id: &str,
        tss_region: &TSSRegion,
    ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
        let mid: i32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

        let mut stmt1: rusqlite::Statement<'_> = pool.prepare(IN_PROMOTER_SQL)?;

        let mapped_rows_1 = stmt1.query_map(
            rusqlite::params![
                mid,
                gene_id,
                location.chr,
                mid,
                tss_region.offset_5p,
                mid,
                tss_region.offset_3p,
            ],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        )?;

        let features_pos =
            mapped_rows_1.filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok());

        let mut stmt2: rusqlite::Statement<'_> = pool.prepare(IN_PROMOTER_SQL)?;

        // negative strand so flip tss region
        let mapped_rows_2 = stmt2.query_map(
            rusqlite::params![
                mid,
                gene_id,
                location.chr,
                mid,
                tss_region.offset_3p,
                mid,
                tss_region.offset_5p,
            ],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        )?;

        let features_neg =
            mapped_rows_2.filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok());

        let features: Vec<GenomicFeature> = features_pos
            .chain(features_neg)
            .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    pub fn get_closest_genes(
        &self,
        location: &dna::Location,
        n: u16,
        level: Level,
    ) -> Result<Vec<GenomicFeature>, Box<dyn Error>> {
        let mid: i32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = self.pool.get()?;

        let mut stmt: rusqlite::Statement<'_> = pool.prepare(CLOSEST_GENE_SQL)?;

        // query_map converts rusqlite into a standard iterator
        let mapped_rows = stmt.query_map(
            rusqlite::params![mid, level as u8, location.chr, mid, n],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        )?;

        // filter map because the query returns an iterator of results
        // and if element is ok, the data is the feature record. Use
        // filter map to keep only the valid records and convert them to
        // actual data by removing the Ok wrapper
        let features: Vec<GenomicFeature> = mapped_rows
            .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
            .collect::<Vec<GenomicFeature>>();

        Ok(features)
    }

    // Returns element
}

// fn mapped_rows_to_features(mapped_rows: MappedRows<'_, impl Fn(&Row<'_>) -> Result<GenomicFeature, Error>>) -> Vec<GenomicFeature> {
//     return mapped_rows
//     .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
//     .collect::<Vec<GenomicFeature>>();
// }

fn row_to_feature(row: &rusqlite::Row<'_>) -> Result<GenomicFeature, rusqlite::Error> {
    let id: u32 = row.get(0)?;
    let chr: String = row.get(1)?;
    let start: i32 = row.get(2)?;
    let end: i32 = row.get(3)?;
    let strand: String = row.get(4)?;
    let gene_id: String = row.get(5)?;
    let gene_symbol: String = row.get(6)?;
    let dist: i32 = row.get(7)?;

    Ok(GenomicFeature {
        id,
        chr,
        start,
        end,
        strand,
        gene_id,
        gene_symbol,
        dist,
    })
}
