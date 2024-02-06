use std::fmt;

use dna::Location;
use r2d2_sqlite::SqliteConnectionManager;

use serde::Serialize;

mod tests;

const WITHIN_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
    FROM genes 
    WHERE level=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) 
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

#[derive(Serialize)]
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

#[derive(Serialize)]
pub struct GenomicFeature {
    pub id: u32,
    pub chr: String,
    pub start: u32,
    pub end: u32,
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
    pub fn new(file: &str) -> Result<Loctogene, String> {
        // let db: Connection = match Connection::open(file) {
        //     Ok(db) => db,
        //     Err(err) => return Err(format!("{}", err)),
        // };

        let manager: SqliteConnectionManager = SqliteConnectionManager::file(file);

        let pool: r2d2::Pool<SqliteConnectionManager> = match r2d2::Pool::builder().build(manager) {
            Ok(pool) => pool,
            Err(err) => return Err(format!("{}", err)),
        };

        Ok(Loctogene { pool })
    }

    pub fn get_genes_within(
        &self,
        location: &Location,
        level: Level,
    ) -> Result<Vec<GenomicFeature>, String> {
        let mid: u32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = match self.pool.get() {
            Ok(pool) => pool,
            Err(err) => return Err(format!("{}", err)),
        };

        let mut stmt: rusqlite::Statement<'_> = match pool.prepare(WITHIN_GENE_SQL) {
            Ok(stmt) => stmt,
            Err(err) => return Err(format!("{}", err)),
        };

        let mapped_rows = match stmt.query_map(
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
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

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
    ) -> Result<Vec<GenomicFeature>, String> {
        let mid: u32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = match self.pool.get() {
            Ok(pool) => pool,
            Err(err) => return Err(format!("{}", err)),
        };

        let mut stmt: rusqlite::Statement<'_> = match pool.prepare(IN_EXON_SQL) {
            Ok(stmt) => stmt,
            Err(err) => return Err(format!("{}", err)),
        };

        let mapped_rows = match stmt.query_map(
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
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

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
    ) -> Result<Vec<GenomicFeature>, String> {
        let mid: u32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = match self.pool.get() {
            Ok(pool) => pool,
            Err(err) => return Err(format!("{}", err)),
        };

        let mut stmt1: rusqlite::Statement<'_> = match pool.prepare(IN_PROMOTER_SQL) {
            Ok(stmt) => stmt,
            Err(err) => return Err(format!("{}", err)),
        };

        let mapped_rows_1 = match stmt1.query_map(
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
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

        let features_pos =
            mapped_rows_1.filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok());

        let mut stmt2: rusqlite::Statement<'_> = match pool.prepare(IN_PROMOTER_SQL) {
            Ok(stmt) => stmt,
            Err(err) => return Err(format!("{}", err)),
        };

        // negative strand so flip tss region
        let mapped_rows_2 = match stmt2.query_map(
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
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

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
    ) -> Result<Vec<GenomicFeature>, String> {
        let mid: u32 = location.mid();

        let pool: r2d2::PooledConnection<SqliteConnectionManager> = match self.pool.get() {
            Ok(pool) => pool,
            Err(err) => return Err(format!("{}", err)),
        };

        let mut stmt: rusqlite::Statement<'_> = match pool.prepare(CLOSEST_GENE_SQL) {
            Ok(stmt) => stmt,
            Err(err) => return Err(format!("{}", err)),
        };

        // query_map converts rusqlite into a standard iterator
        let mapped_rows = match stmt.query_map(
            rusqlite::params![mid, level as u8, location.chr, mid, n],
            |row: &rusqlite::Row<'_>| row_to_feature(row),
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

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
    let id: u32 = match row.get(0) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let chr: String = match row.get(1) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let start: u32 = match row.get(2) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let end: u32 = match row.get(3) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let strand: String = match row.get(4) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let gene_id: String = match row.get(5) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let gene_symbol: String = match row.get(6) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

    let dist: i32 = match row.get(7) {
        Ok(v) => v,
        Err(err) => return Err(err),
    };

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
