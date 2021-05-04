use seine::salmon::{EqClassCollection, EqClassView, FromPathExt, QuantEntry};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::path::Path;

use clap::{crate_authors, crate_version, value_t};
use clap::{App, Arg};
use serde::Serialize;
use slog::{info, o, Drain};
use std::fs;

fn main() -> Result<(), serde_yaml::Error> {
    let crate_authors = crate_authors!("\n");
    let version = crate_version!();

    // Logging
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator)
        .use_custom_timestamp(|out: &mut dyn std::io::Write| {
            write!(out, "{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")).unwrap();
            Ok(())
        })
        .build()
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    let log = slog::Logger::root(drain, o!());

    let eval_app = App::new("eval")
        .version(version)
        .author(crate_authors)
        .arg(
            Arg::from_usage("--train_ecs=<train_ecs> 'Training set equivalence class counts'")
                .required(true),
        )
        .arg(Arg::from_usage("--quants=<train_ecs> 'Training set quants'").required(true))
        .arg(
            Arg::from_usage("--test_ecs=<test_salmon_dir> 'Test set equivalence class counts'")
                .required(true),
        )
        .arg(
            Arg::from_usage("--smoothing=<smoothing> 'TPM based smoothing'")
                .required(false)
                // .default_value("1"),
                .default_value("1e-8"),
        )
        .arg(Arg::from_usage(
            "-t, --smooth_per_tx 'use reads per tx smoothing'",
        ))
        .arg(Arg::from_usage("-o --output=<output> 'output results path'").required(true));

    let opts = App::new("perplexity")
        .version(version)
        .author(crate_authors)
        .subcommand(eval_app)
        .get_matches();

    if let Some(opts) = opts.subcommand_matches("eval") {
        let output_file = opts.value_of("output").unwrap().to_string();
        let train_ecs = opts.value_of("train_ecs").unwrap().to_string();
        let train_quants = opts.value_of("quants").unwrap().to_string();
        let test_ecs = opts.value_of("test_ecs").unwrap().to_string();
        let alpha_smooth = value_t!(opts.value_of("smoothing"), f64).unwrap();

        let smooth_per_tx = opts.is_present("smooth_per_tx");

        let output_file = Path::new(&output_file);

        let mut file_name = String::from(output_file.file_stem().unwrap().to_str().unwrap());
        file_name.push_str("_auxinfo.csv");
        let auxinfo_file = output_file.parent().unwrap().join(file_name);

        if smooth_per_tx {
            info!(log, "Smoothing reads per transcript");
        } else {
            info!(log, "Smoothing per TPM";)
        }
        let result = perplexity_from_files(
            &train_quants,
            &train_ecs,
            &test_ecs,
            alpha_smooth,
            smooth_per_tx,
            &log,
        );
        let summary = perplexity_summary(&result);

        info!(log, "{:#?}", &summary);

        info!(log, "Writing auxinfo to: {:?}", &auxinfo_file);
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b',')
            .from_path(auxinfo_file)
            .unwrap();

        for record in result.results {
            wtr.serialize(record).unwrap();
        }

        info!(log, "Writing perplexity result to: {:?}", output_file);
        let s = serde_yaml::to_string(&summary)?;
        fs::write(output_file, s).expect("Unable to write file");
    }
    Ok(())
}

fn perplexity_from_files<P: AsRef<Path>>(
    quants: &P,
    train_ecs: &P,
    test_ecs: &P,
    alpha_smooth: f64,
    smooth_per_tx: bool,
    log: &slog::Logger,
) -> EqClassCollectionEval {
    info!(log, "Loading quants from {:?}", &quants.as_ref());
    let quant_map = HashMap::<String, QuantEntry>::from_path(&quants).unwrap();

    info!(log, "Loading quantified ECs from {:?}", &train_ecs.as_ref());
    let tr_ecs = EqClassCollection::from_path(&train_ecs).unwrap();

    info!(log, "Loading held-out ECs from {:?}", &test_ecs.as_ref());
    let te_ecs = EqClassCollection::from_path(&test_ecs).unwrap();

    let ec_eval_results = perplexity(
        te_ecs,
        &quant_map,
        &tr_ecs,
        alpha_smooth,
        smooth_per_tx,
        &log,
    );
    ec_eval_results
}

pub fn mapped_txs(ecs: &EqClassCollection) -> HashSet<String> {
    let mut set = HashSet::new();
    for ec in ecs.classes.iter() {
        for label in ec.labels {
            set.insert(ecs.targets[*label].clone());
        }
    }
    set
}

pub fn perplexity(
    ecs: EqClassCollection,
    quant_map: &HashMap<String, QuantEntry>,
    train_ecs: &EqClassCollection,
    alpha_smooth: f64,
    smooth_per_tx: bool,
    log: &slog::Logger,
) -> EqClassCollectionEval {
    let train_txs = mapped_txs(train_ecs);

    let smoothed_quant_map = if smooth_per_tx {
        per_tx_smoothing(quant_map, alpha_smooth)
    } else {
        per_TPM_smoothing(quant_map, alpha_smooth)
    };

    let mut results = Vec::new();
    // For each EC
    for ec in ecs.classes.iter() {
        let EqClassView {
            labels,
            weights,
            count,
        } = ec;

        let targets = &ecs.targets;
        let label_names: Vec<String> = labels.iter().map(|l| targets[*l].clone()).collect();
        let mut ec_perp = 0.;

        let mut ec_perp_smooth = 0.;

        // Sum out conditional probs per tx
        for (l, w) in labels.iter().zip(weights.iter()) {
            let t = &ecs.targets[*l];

            // not-smoothed probability for impossible read analysis
            let eta = quant_map[t].tpm / 1e6;
            ec_perp += w * eta;

            // smoothed probability
            let eta = smoothed_quant_map[t].tpm / 1e6;
            ec_perp_smooth += w * eta;
        }

        let class = held_out_ec_class(ec_perp, &label_names, &train_txs);

        let log_l = if ec_perp > 0. {
            Some(ec_perp.ln())
        } else {
            None
        };

        let log_l_smooth = ec_perp_smooth.ln();

        let count = count as usize;
        let eval_info = EqClassEvalInfo {
            log_l,
            log_l_smooth,
            count,
            class,
        };
        results.push(eval_info);
    }

    EqClassCollectionEval { ecs, results }
}

pub fn held_out_ec_class(
    prob: f64,
    labels: &[String],
    train_txs: &HashSet<String>,
) -> EqClassEvalClass {
    // Ambig or unique
    if labels.len() == 1 {
        let tx = &labels[0];
        if prob > 0. {
            assert!(train_txs.contains(tx));
            EqClassEvalClass::PossibleUnique
        } else {
            if train_txs.contains(tx) {
                EqClassEvalClass::ImpossibleUnique
            } else {
                EqClassEvalClass::ImpossibleUniqueNovel
            }
        }
    } else if labels.len() > 1 {
        if prob > 0. {
            if labels.iter().all(|t| train_txs.contains(t)) {
                EqClassEvalClass::PossibleAmbig
            } else {
                EqClassEvalClass::PossibleAmbigNovel
            }
        } else {
            // contains novel, only novel, or no novel
            if labels.iter().any(|t| train_txs.contains(t)) {
                if labels.iter().all(|t| train_txs.contains(t)) {
                    EqClassEvalClass::ImpossibleAmbig
                } else {
                    EqClassEvalClass::ImpossibleAmbigContainsNovel
                }
            } else {
                EqClassEvalClass::ImpossibleAmbigOnlyNovel
            }
        }
    } else {
        panic!("Zero length EC label");
    }
}

pub fn per_tx_smoothing(
    quant_map: &HashMap<String, QuantEntry>,
    alpha: f64,
) -> HashMap<String, QuantEntry> {
    let mut smoothed_quants = HashMap::new();

    // compute the new denominator
    let mut denom = 0.0;
    for (_, v) in quant_map.iter() {
        denom += (v.num_reads + alpha) / v.efflen;
    }

    for (k, v) in quant_map.iter() {
        let num_reads = v.num_reads + alpha;
        let tpm = 1e6 * num_reads / v.efflen / denom;
        let smoothed_entry = QuantEntry {
            tpm,
            num_reads,
            efflen: v.efflen,
            len: v.len,
        };
        smoothed_quants.insert(k.clone(), smoothed_entry);
    }
    smoothed_quants
}

pub fn per_TPM_smoothing(
    quant_map: &HashMap<String, QuantEntry>,
    alpha: f64,
) -> HashMap<String, QuantEntry> {
    let mut smoothed_quants = HashMap::new();

    let m = quant_map.len();

    let denom = 1. + (m as f64 * alpha);

    // reconstituted counts
    let mut total_reads = 0.;
    for (_, v) in quant_map.iter() {
        total_reads = total_reads + v.num_reads;
    }

    for (k, v) in quant_map.iter() {
        let eta = ((v.tpm / 1e6) + alpha) / denom;
        let tpm = eta * 1e6;
        let num_reads = total_reads * eta;

        let smoothed_entry = QuantEntry {
            tpm,
            num_reads,
            efflen: v.efflen,
            len: v.len,
        };
        smoothed_quants.insert(k.clone(), smoothed_entry);
    }

    smoothed_quants
}

pub fn perplexity_summary(results: &EqClassCollectionEval) -> ResultSummary {
    let n_ecs = results.ecs.neq;
    let n_reads: usize = results.results.iter().map(|i| i.count).sum();

    let n_possible_ecs: usize = results
        .results
        .iter()
        .filter(|i| is_possible(i.class))
        .count();
    let n_possible_reads: usize = results
        .results
        .iter()
        .filter(|i| is_possible(i.class))
        .map(|i| i.count)
        .sum();

    let n_impossible_ecs: usize = results
        .results
        .iter()
        .filter(|i| !is_possible(i.class))
        .count();
    let n_impossible_reads: usize = results
        .results
        .iter()
        .filter(|i| !is_possible(i.class))
        .map(|i| i.count)
        .sum();

    let n_discarded_ecs: usize = results
        .results
        .iter()
        .filter(|i| is_discarded(i.class))
        .count();
    let n_discarded_reads: usize = results
        .results
        .iter()
        .filter(|i| is_discarded(i.class))
        .map(|i| i.count)
        .sum();

    // perplexity with smoothing.
    let mut tot_reads = 0.;
    let mut tot_perp = 0.;
    for ec_info in results.results.iter() {
        let class = ec_info.class;
        let count = ec_info.count as f64;
        if !is_discarded(class) {
            tot_reads += count;
            tot_perp += count * ec_info.log_l_smooth;
        }
    }
    let smoothed_perplexity = -tot_perp / tot_reads;
    let smoothed_perplexity = smoothed_perplexity.exp2();

    let class_ec_freq: Vec<EqClassEvalClass> = results.results.iter().map(|i| i.class).collect();
    let class_ec_freq = to_freq_table(&class_ec_freq);
    let mut class_read_freq = HashMap::new();
    let mut class_read_count_freq = HashMap::new();

    for i in results.results.iter() {
        let count = class_read_freq.entry(i.class).or_insert(0);
        *count += i.count;

        if !is_possible(i.class) {
            let freq_table = class_read_count_freq
                .entry(i.class)
                .or_insert(HashMap::new());
            let count = freq_table.entry(i.count).or_insert(0);
            *count += 1;
        }
    }

    ResultSummary {
        smoothed_perplexity,
        n_ecs,
        n_reads,
        n_possible_ecs,
        n_possible_reads,
        n_impossible_ecs,
        n_impossible_reads,
        n_discarded_ecs,
        n_discarded_reads,

        class_ec_freq,
        class_read_freq,
        class_read_count_freq,
    }
}

pub fn to_freq_table<T>(arr: &[T]) -> HashMap<T, usize>
where
    T: Hash + Eq + Clone,
{
    let mut freqtable = HashMap::new();
    for item in arr.iter() {
        let counter = freqtable.entry(item.clone()).or_insert(0);
        *counter += 1;
    }
    freqtable
}

#[derive(Debug, Serialize, Default)]
pub struct ResultSummary {
    smoothed_perplexity: f64,
    n_ecs: usize,
    n_reads: usize,
    n_possible_ecs: usize,
    n_possible_reads: usize,
    n_impossible_ecs: usize,
    n_impossible_reads: usize,
    n_discarded_ecs: usize,
    n_discarded_reads: usize,

    class_ec_freq: HashMap<EqClassEvalClass, usize>,
    class_read_freq: HashMap<EqClassEvalClass, usize>,
    class_read_count_freq: HashMap<EqClassEvalClass, HashMap<usize, usize>>,
}

#[derive(Debug, Serialize, PartialEq, Eq, Hash, Copy, Clone)]
pub enum EqClassEvalClass {
    PossibleUnique, // Uniquely mapped to tx observed at training time
    //PossibleUniqueNovel     <-- not possible!
    PossibleAmbig,      // Ambiguously mapped to txs observed at trainging time
    PossibleAmbigNovel, // Ambigously mapped at least one novel tx
    //PossibleAmbigOnlyNovel, <-- not possible!
    ImpossibleUnique, // uniquely mapped read for wich tx is ambig only in training
    ImpossibleUniqueNovel, // uniquely mapped to novel tx
    ImpossibleAmbig,  // ambig mapped reads to new ec pattern containing no novel txs
    ImpossibleAmbigContainsNovel, // // Ambiguously mapped to at least one novel and one observed tx
    ImpossibleAmbigOnlyNovel, // Ambiguously mapped only to novel txs
    PANIC,
}

pub fn is_possible(class: EqClassEvalClass) -> bool {
    matches!(
        class,
        EqClassEvalClass::PossibleUnique
            | EqClassEvalClass::PossibleAmbig
            | EqClassEvalClass::PossibleAmbigNovel
    )
}

pub fn is_discarded(class: EqClassEvalClass) -> bool {
    matches!(
        class,
        EqClassEvalClass::ImpossibleUniqueNovel | EqClassEvalClass::ImpossibleAmbigOnlyNovel
    )
}

#[derive(Debug, Serialize)]
pub struct EqClassEvalInfo {
    #[serde(rename = "Log_L")]
    log_l: Option<f64>, //None if probability == 0

    #[serde(rename = "Log_L_smoothed")]
    log_l_smooth: f64,

    #[serde(rename = "Count")]
    count: usize,

    #[serde(rename = "Class")]
    class: EqClassEvalClass,
}

impl Default for EqClassEvalClass {
    fn default() -> Self {
        EqClassEvalClass::PANIC
    }
}

#[derive(Debug)]
pub struct EqClassCollectionEval {
    ecs: EqClassCollection,
    results: Vec<EqClassEvalInfo>,
}
