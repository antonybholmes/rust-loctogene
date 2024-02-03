use rusqlite::Connection;
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

#[derive(Serialize)]
pub struct FeatureRecord {
    pub id: u32,
    pub chr: String,
    pub start: u32,
    pub end: u32,
    pub strand: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub dist: i32,
}

#[derive(Serialize)]
pub struct Features {
    pub location: String,
    pub level: String,
    pub features: Vec<FeatureRecord>,
}

//const NO_FEATURES: [Features; 0] = [] .to_vec();

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

fn get_level(level: u32) -> String {
    match level {
        2 => "transcript".to_string(),
        3 => "exon".to_owned(),
        _ => "gene".to_owned(),
    }
}

pub struct Loctogene {
    db: Connection,
}

impl Loctogene {
    pub fn new(file: &str) -> Result<Loctogene, String> {
        let db: Connection = Connection::open(file).unwrap();

        return Ok(Loctogene { db });
    }

    pub fn get_genes_within(
        &self,
        location: &dna::Location,
        level: u32,
    ) -> Result<Features, String> {
        let mid: u32 = (location.start + location.end) / 2;

        let mut stmt = match self.db.prepare(WITHIN_GENE_SQL) {
            Ok(stmt) => stmt,
            Err(err) => panic!("{}", err),
        };

        let mapped_rows = match stmt.query_map(
            rusqlite::params![
                mid,
                level,
                location.chr,
                location.start,
                location.start,
                location.end,
                location.end
            ],
            |row| {
                Ok(FeatureRecord {
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
            Err(err) => panic!("{}", err),
        };

        // let mut id:i64 = 0;
        // while let Some(row) = rows.next().expect("while row failed") {
        //     id=row.get(1).expect("get row failed");
        // }

        let records: Vec<FeatureRecord> = mapped_rows
            .filter_map(Result::ok)
            .collect::<Vec<FeatureRecord>>();

        let ret: Features = Features {
            location: format!("{}", location),
            level: get_level(level),
            features: records,
        };

        return Ok(ret);
    }

    pub fn get_closest_genes(
        &self,
        location: &dna::Location,
        n: u16,
        level: u32,
    ) -> Result<Features, String> {
        let mid: u32 = (location.start + location.end) / 2;

        let mut stmt = match self.db.prepare(CLOSEST_GENE_SQL) {
            Ok(stmt) => stmt,
            Err(err) => panic!("{}", err),
        };

        let mapped_rows =
            match stmt.query_map(rusqlite::params![mid, level, location.chr, mid, n], |row| {
                Ok(FeatureRecord {
                    id: row.get(0).expect("col 0"),
                    chr: row.get(1).expect("col 1"),
                    start: row.get(2).expect("col 2"),
                    end: row.get(3).expect("col 3"),
                    strand: row.get(4).expect("col 4"),
                    gene_id: row.get(5).expect("col 5"),
                    gene_symbol: row.get(6).expect("col 6"),
                    dist: row.get(7).expect("col 7"),
                })
            }) {
                Ok(rows) => rows,
                Err(err) => panic!("{}", err),
            };

        let records: Vec<FeatureRecord> = mapped_rows
            .filter_map(Result::ok)
            .collect::<Vec<FeatureRecord>>();

        let ret: Features = Features {
            location: format!("{}", location),
            level: get_level(level),
            features: records,
        };

        return Ok(ret);
    }
}
