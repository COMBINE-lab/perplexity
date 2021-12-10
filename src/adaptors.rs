use clap::{value_t, ArgMatches};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};

use seine::salmon::QuantEntry;

/*******************************************************************************/
/*                         Quants                                              */
/*******************************************************************************/

#[derive(Debug, Deserialize, Serialize)]
struct EXpressRecord {
    pub bundle_id: u32,
    pub target_id: String,
    pub length: u32,
    pub eff_length: f64,
    pub tot_counts: u32,
    pub uniq_counts: u32,
    pub est_counts: f64,
    pub eff_counts: f64,
    pub ambig_distr_alpha: f64,
    pub ambig_distr_beta: f64,
    pub fpkm: f64,
    pub fpkm_conf_low: f64,
    pub fpkm_conf_high: f64,
    pub solvable: char, //TODO implement deserializer for "T" and "F"
    pub tpm: f64,
}

/*******************************************************************************/
/*                         Extension Traits                                    */
/*******************************************************************************/

// For express
trait FromEXpressPathExt {
    fn from_express_path<P: AsRef<Path>>(p: P) -> Result<Self, csv::Error>
    where
        Self: Sized;
}

impl FromEXpressPathExt for HashMap<String, QuantEntry> {
    fn from_express_path<P: AsRef<Path>>(p: P) -> Result<HashMap<String, QuantEntry>, csv::Error> {
        // TODO error handling for csv::Error and io::Error
        let file = File::open(p);
        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file.unwrap());

        let mut quant_map = HashMap::new();

        for quant_record in rdr.deserialize() {
            let express_record: EXpressRecord = quant_record?;
            let efflen = express_record.eff_length;

            // A necssary evil? eXpress does not output efflens of unxpressed txps
            let efflen = if efflen > 0. {
                efflen
            } else {
                express_record.length as f64
            };
            let quant_entry = QuantEntry {
                len: express_record.length,
                efflen: efflen,
                tpm: express_record.tpm,
                num_reads: express_record.eff_counts,
            };

            quant_map.insert(express_record.target_id.to_string(), quant_entry);
        }

        Ok(quant_map)
    }
}

pub fn load_express_quants<P: AsRef<Path>>(
    p: P,
) -> Result<HashMap<String, QuantEntry>, csv::Error> {
    HashMap::<String, QuantEntry>::from_express_path(p)
}

pub fn load_express_txps_w_mapped_reads<P: AsRef<Path>>(
    p: P,
) -> Result<HashSet<String>, csv::Error> {
    // TODO error handling for csv::Error and io::Error
    let file = File::open(p);
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file.unwrap());

    let mut txps_w_mapped_reads = HashSet::new();

    for quant_record in rdr.deserialize() {
        let express_record: EXpressRecord = quant_record?;
        if express_record.tot_counts > 0 {
            txps_w_mapped_reads.insert(express_record.target_id.clone());
        }
    }

    Ok(txps_w_mapped_reads)
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ExpressAppInfo {
    // Arguments
    pub mode: String,
    pub output_path: PathBuf,
    pub quant_path: PathBuf,
    pub val_ecs_path: PathBuf,
    pub smoothing_strategy: String, // TODO maybe an enum?
    pub smoothing_param: Option<f64>,

    // Extra stuff
    pub auxinfo_path: PathBuf,
    pub appinfo_path: PathBuf,
}

impl ExpressAppInfo {
    pub fn from_matches(matches: &ArgMatches) -> Self {
        let output_path = matches.value_of("output").unwrap().to_string();

        let quant_path = matches.value_of("quants").unwrap().to_string();
        let quant_path = PathBuf::from(&quant_path);
        let val_ecs_path = matches.value_of("test_ecs").unwrap().to_string();
        let val_ecs_path = PathBuf::from(&val_ecs_path);

        let smoothing_strategy = matches.value_of("smoothing_strategy").unwrap().to_string();

        let smoothing_param = match smoothing_strategy.as_str() {
            "TX" | "TPM" => {
                let param = value_t!(matches.value_of("smoothing"), f64).unwrap();
                Some(param)
            }
            _ => None,
        };

        let output_path = PathBuf::from(&output_path);

        let prefix = String::from(output_path.file_stem().unwrap().to_str().unwrap());

        let mut auxinfo_filename = prefix.clone();
        auxinfo_filename.push_str("_auxinfo.csv");
        let auxinfo_path = output_path.parent().unwrap().join(auxinfo_filename);

        let mut appinfo_filename = prefix.clone();
        appinfo_filename.push_str("_appinfo.yml");
        let appinfo_path = output_path.parent().unwrap().join(appinfo_filename);

        Self {
            mode: String::from("eXpress"),
            output_path,
            quant_path,
            val_ecs_path,
            smoothing_strategy,
            smoothing_param,

            // Extra stuff
            auxinfo_path,
            appinfo_path,
        }
    }
}
