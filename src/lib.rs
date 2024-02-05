use std::fmt;

use dna::Location;
use r2d2_sqlite::SqliteConnectionManager;
use serde::Serialize;

mod tests;

const WITHIN_GENE_SQL: &str = r#"SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? 
    FROM genes 
    WHERE level=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) 
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
            |row: &rusqlite::Row<'_>| {
                Ok(GenomicFeature {
                    id: row.get(0).expect("col 0"),
                    chr: row.get(1).expect("col 1"),
                    start: row.get(2).expect("col 2"),
                    end: row.get(3).expect("col 3"),
                    strand: row.get(4).expect("col 4"),
                    gene_id: row.get(5).expect("col 5"),
                    gene_symbol: row.get(6).expect("col 6"),
                    dist: row.get(7).expect("col 7"),
                })
            },
        ) {
            Ok(rows) => rows,
            Err(err) => return Err(format!("{}", err)),
        };

        let features: Vec<GenomicFeature> = mapped_rows
            .filter_map(|x: Result<GenomicFeature, rusqlite::Error>| x.ok())
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

        println!("{}", level as u8);

        // query_map converts rusqlite into a standard iterator
        let mapped_rows = match stmt.query_map(
            rusqlite::params![mid, level as u8, location.chr, mid, n],
            |row: &rusqlite::Row<'_>| {
                Ok(GenomicFeature {
                    id: row.get(0).expect("col 0"),
                    chr: row.get(1).expect("col 1"),
                    start: row.get(2).expect("col 2"),
                    end: row.get(3).expect("col 3"),
                    strand: row.get(4).expect("col 4"),
                    gene_id: row.get(5).expect("col 5"),
                    gene_symbol: row.get(6).expect("col 6"),
                    dist: row.get(7).expect("col 7"),
                })
            },
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
}
