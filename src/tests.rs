#[cfg(test)]
use dna::Location;
#[cfg(test)]
use serde_json::json;

#[cfg(test)]
use crate::Loctogene;



#[test]
fn test_within() {
 
    let loc: Location = match Location::parse("chr3:187721370-187733550") {
        Ok(loc)=>loc,
        Err(err)=>panic!("{}", err)
    };

    let genesdb: Loctogene = match Loctogene::new("../docker-rust-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records:Vec<crate::GenomicFeature>  =  match genesdb.get_genes_within(&loc, crate::Level::Gene) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js: serde_json::Value = json!(records);

    println!("{}", js);

}

#[test]
fn test_closest() {
    let loc: Location = match Location::parse("chr3:187721370-187733550") {
        Ok(loc)=>loc,
        Err(err)=>panic!("{}", err)
    };

    let genesdb: Loctogene = match Loctogene::new("../docker-rust-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records:Vec<crate::GenomicFeature>  =  match genesdb.get_closest_genes(&loc, 10, crate::Level::Gene) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js: serde_json::Value = json!(records);

    println!("{}", js);

}