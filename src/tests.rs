#[cfg(test)]
use dna::Location;
#[cfg(test)]
use serde_json::json;

#[cfg(test)]
use crate::Loctogene;



#[test]
fn test_within() {
    let loc: Location= Location::parse("chr3:187721370-187733550");

    let genesdb: Loctogene = match Loctogene::new("/home/antony/development/go/docker-go-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records =  match genesdb.get_genes_within(&loc, 1) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js = json!(records);

    println!("{}", js);

}

#[test]
fn test_closest() {
    let loc: Location= Location::parse("chr3:187721377-187745725");

    let genesdb: Loctogene = match Loctogene::new("/home/antony/development/go/docker-go-edb-api/data/loctogene/grch38.db") {
        Ok(db)=>db,
        Err(err)=>panic!("{}", err)
    };

    let records =  match genesdb.get_closest_genes(&loc, 10, 1) {
        Ok(records)=>records,
        Err(err)=>panic!("{}", err)
    };

    let js = json!(records);

    println!("{}", js);

}