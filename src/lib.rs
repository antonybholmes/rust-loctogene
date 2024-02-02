use rusqlite::Connection;
use serde::Serialize;

mod tests;

const WITHIN_GENE_SQL: &str = "SELECT id, chr, start, end, strand, gene_id, gene_symbol, start - ? FROM genes WHERE level=? AND chr=? AND ((start <= ? AND end >= ?) OR (start <= ? AND end >= ?)) ORDER BY start ASC";

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

        let mut rows: rusqlite::Rows<'_> = match stmt.query(rusqlite::params![
            mid,
            level,
            location.chr,
            location.start,
            location.start,
            location.end,
            location.end
        ]) {
            Ok(rows) => rows,
            Err(err) => panic!("{}", err),
        };

        // let mut id:i64 = 0;
        // while let Some(row) = rows.next().expect("while row failed") {
        //     id=row.get(1).expect("get row failed");
        // }

        let ret: Features = match rows_to_records(&location, &mut rows, level) {
            Ok(features) => features,
            Err(err) => panic!("{}", err),
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

        let mut rows: rusqlite::Rows<'_> =
            match stmt.query(rusqlite::params![mid, level, location.chr, mid, n]) {
                Ok(rows) => rows,
                Err(err) => panic!("{}", err),
            };

        let ret: Features = match rows_to_records(&location, &mut rows, level) {
            Ok(features) => features,
            Err(err) => panic!("{}", err),
        };

        return Ok(ret);
    }
}

 
fn rows_to_records(
    location: &dna::Location,
    rows: &mut rusqlite::Rows<'_>,
    level: u32,
) -> Result<Features, String> {
    // let id: u32;
    // let chr: String;
    // let start: u32;
    // let end: u32;
    // let strand: String;
    // let gene_id: String;
    // let gene_symbol: String;
    // let d: u32;

    let t: String = match level {
        2 => "transcript".to_string(),
        3 => "exon".to_owned(),
        _ => "gene".to_owned(),
    };

    let mut records: Vec<FeatureRecord> = Vec::new();

    while let Some(row) = rows.next().expect("while row failed") {
        // if strand == "-" {
        //     t := start
        //     start = end
        //     end = t
        // }

        let record = FeatureRecord {
            id: row.get(0).expect("0 row failed"),
            chr: row.get(1).expect("1 row failed"),
            start: row.get(2).expect("2 row failed"),
            end: row.get(3).expect("3 row failed"),
            strand: row.get(4).expect("4 row failed"),
            gene_id: row.get(5).expect("5 row failed"),
            gene_symbol: row.get(6).expect("6 row failed"),
            dist: row.get(7).expect("7 row failed"),
        };

        records.push(record);
    }

    let ret: Features = Features {
        location: format!("{}", location),
        level: t.to_string(),
        features: records,
    };

    return Ok(ret);
}
