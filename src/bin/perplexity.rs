use seine::salmon::{EqClassCollection, EqClassView, FromPathExt, QuantEntry};
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::path::PathBuf;

use clap::{crate_authors, crate_version, value_t};
use clap::{App, Arg, ArgMatches};
use serde::Serialize;
use slog::{info, o, Drain};
use std::fs;

use perplexity::adaptors::{load_express_quants, load_express_txps_w_mapped_reads, ExpressAppInfo};
use perplexity::smoothing::{PerTPMSmoother, PerTxpSmoother, SGTSmoother, SmoothingStrategy};

/*******************************************************************************/
// App
/*******************************************************************************/

fn app() -> App<'static, 'static> {
    let crate_authors = crate_authors!("\n");
    let version = crate_version!();

    let eval_app = App::new("eval")
        .version(version)
        .author(crate_authors)
        .arg(
            Arg::from_usage("--train_ecs=<train_ecs> 'Training set equivalence class file'")
                .required(true),
        )
        .arg(Arg::from_usage("--quants=<train_ecs> 'Training set quant.sf'").required(true))
        .arg(
            Arg::from_usage("--test_ecs=<test_salmon_dir> 'Test set equivalence class file'")
                .required(true),
        )
        .arg(
            Arg::from_usage("--smoothing=<smoothing> 'Smoothing parameter for various strategies'")
                .required(false)
                .default_value("1e-8"),
        )
        .arg(
            Arg::from_usage("-m, --smoothing_strategy=<mode> 'One of {TX, TPM, LGT}'")
                .required(false)
                .default_value("TPM"),
        )
        .arg(Arg::from_usage("-o --output=<output> 'Write results to output'").required(true));

    let xprs_app = App::new("eval-xprs")
        .version(version)
        .author(crate_authors)
        .arg(Arg::from_usage("--quants=<train_ecs> 'Training set quant.sf'").required(true))
        .arg(
            Arg::from_usage("--test_ecs=<test_salmon_dir> 'Test set equivalence class file'")
                .required(true),
        )
        .arg(
            Arg::from_usage("--smoothing=<smoothing> 'Smoothing parameter for various strategies'")
                .required(false)
                .default_value("1e-8"),
        )
        .arg(
            Arg::from_usage("-m, --smoothing_strategy=<mode> 'One of {TX, TPM, LGT}'")
                .required(false)
                .default_value("TPM"),
        )
        .arg(Arg::from_usage("-o --output=<output> 'Write results to output'").required(true));

    // return the app
    App::new("perplexity")
        .version(version)
        .author(crate_authors)
        .subcommand(eval_app)
        .subcommand(xprs_app)
}

#[derive(Debug, Serialize)]
struct AppInfo {
    // Arguments
    output_path: PathBuf,
    quant_path: PathBuf,
    quant_ecs_path: PathBuf,
    val_ecs_path: PathBuf,
    smoothing_strategy: String, // TODO maybe an enum?
    smoothing_param: Option<f64>,

    // Extra stuff
    auxinfo_path: PathBuf,
    appinfo_path: PathBuf,
}

impl AppInfo {
    fn from_matches(matches: &ArgMatches) -> Self {
        let output_path = matches.value_of("output").unwrap().to_string();

        let quant_path = matches.value_of("quants").unwrap().to_string();
        let quant_path = PathBuf::from(&quant_path);
        let quant_ecs_path = matches.value_of("train_ecs").unwrap().to_string();
        let quant_ecs_path = PathBuf::from(&quant_ecs_path);
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
            output_path,
            quant_path,
            quant_ecs_path,
            val_ecs_path,
            smoothing_strategy,
            smoothing_param,

            // Extra stuff
            auxinfo_path,
            appinfo_path,
        }
    }
}

/*******************************************************************************/
// main
/*******************************************************************************/

fn main() -> Result<(), serde_yaml::Error> {
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

    // App and parse args
    let app = app();
    let matches = app.get_matches();

    if let Some(matches) = matches.subcommand_matches("eval") {
        let app_info = AppInfo::from_matches(&matches);
        let result = perplexity_from_app_info(&app_info, &log);
        let summary = perplexity_summary(&result);

        info!(log, "Writing auxinfo to: {:?}", &app_info.auxinfo_path);
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b',')
            .from_path(&app_info.auxinfo_path)
            .unwrap();

        for record in result.results {
            wtr.serialize(record).unwrap();
        }

        info!(log, "Writing AppInfo to: {:?}", &app_info.appinfo_path);
        let app_info_yaml = serde_yaml::to_string(&app_info)?;
        fs::write(&app_info.appinfo_path, app_info_yaml).expect("Unable to write file");
        info!(
            log,
            "Writing perplexity result to: {:?}", &app_info.output_path
        );
        let s = serde_yaml::to_string(&summary)?;
        fs::write(&app_info.output_path, s).expect("Unable to write file");

        info!(log, "# Total validation reads: {:?}", &summary.n_reads);
        info!(log, "# Impossible reads: {:?}", &summary.n_impossible_reads);
        info!(log, "# Discarded reads: {:?}", &summary.n_discarded_reads);
        info!(log, "Perplexity: {:?}", &summary.smoothed_perplexity);
    }

    if let Some(matches) = matches.subcommand_matches("eval-xprs") {
        let app_info = ExpressAppInfo::from_matches(&matches);
        let result = perplexity_from_xprs_app_info(&app_info, &log);
        let summary = perplexity_summary(&result);

        info!(log, "Writing auxinfo to: {:?}", &app_info.auxinfo_path);
        let mut wtr = csv::WriterBuilder::new()
            .delimiter(b',')
            .from_path(&app_info.auxinfo_path)
            .unwrap();

        for record in result.results {
            wtr.serialize(record).unwrap();
        }

        info!(log, "Writing AppInfo to: {:?}", &app_info.appinfo_path);
        let app_info_yaml = serde_yaml::to_string(&app_info)?;
        fs::write(&app_info.appinfo_path, app_info_yaml).expect("Unable to write file");
        info!(
            log,
            "Writing perplexity result to: {:?}", &app_info.output_path
        );
        let s = serde_yaml::to_string(&summary)?;
        fs::write(&app_info.output_path, s).expect("Unable to write file");

        info!(log, "# Total validation reads: {:?}", &summary.n_reads);
        info!(log, "# Impossible reads: {:?}", &summary.n_impossible_reads);
        info!(log, "# Discarded reads: {:?}", &summary.n_discarded_reads);
        info!(log, "Perplexity: {:?}", &summary.smoothed_perplexity);
    }

    Ok(())
}

fn perplexity_from_app_info(app_info: &AppInfo, log: &slog::Logger) -> EqClassCollectionEval {
    let smoothing_strategy = &app_info.smoothing_strategy;
    info!(log, "Smoothing strategy {:?}", &smoothing_strategy);

    let smoothing_strategy = {
        match app_info.smoothing_strategy.as_str() {
            "TX" => {
                let smoothing_param = app_info.smoothing_param.unwrap();
                info!(log, "Smoothing parameter {:?}", &smoothing_param);
                Some(Box::new(PerTxpSmoother::new(smoothing_param)) as Box<dyn SmoothingStrategy>)
            }
            "TPM" => {
                let smoothing_param = app_info.smoothing_param.unwrap();
                info!(log, "Smoothing parameter {:?}", &smoothing_param);
                Some(Box::new(PerTPMSmoother::new(smoothing_param)) as Box<dyn SmoothingStrategy>)
            }
            "LGT" => Some(Box::new(SGTSmoother::new()) as Box<dyn SmoothingStrategy>),
            _ => None,
        }
    };

    let smoothing_strategy = smoothing_strategy.unwrap();

    info!(log, "Loading quants from {:?}", &app_info.quant_path);
    let quant_map = HashMap::<String, QuantEntry>::from_path(&app_info.quant_path).unwrap();

    info!(
        log,
        "Loading quantified ECs from {:?}", &app_info.quant_ecs_path
    );
    let tr_ecs = EqClassCollection::from_path(&app_info.quant_ecs_path).unwrap();
    let txs_w_mapped_reads = mapped_txs(&tr_ecs);

    info!(
        log,
        "Loading held-out ECs from {:?}", &app_info.val_ecs_path
    );
    let te_ecs = EqClassCollection::from_path(&app_info.val_ecs_path).unwrap();

    let ec_eval_results = perplexity(
        te_ecs,
        &quant_map,
        &txs_w_mapped_reads,
        smoothing_strategy.as_ref(),
        &log,
    );
    ec_eval_results
}

fn perplexity_from_xprs_app_info(
    app_info: &ExpressAppInfo,
    log: &slog::Logger,
) -> EqClassCollectionEval {
    let smoothing_strategy = &app_info.smoothing_strategy;
    info!(log, "Smoothing strategy {:?}", &smoothing_strategy);

    let smoothing_strategy = {
        match app_info.smoothing_strategy.as_str() {
            "TX" => {
                let smoothing_param = app_info.smoothing_param.unwrap();
                info!(log, "Smoothing parameter {:?}", &smoothing_param);
                Some(Box::new(PerTxpSmoother::new(smoothing_param)) as Box<dyn SmoothingStrategy>)
            }
            "TPM" => {
                let smoothing_param = app_info.smoothing_param.unwrap();
                info!(log, "Smoothing parameter {:?}", &smoothing_param);
                Some(Box::new(PerTPMSmoother::new(smoothing_param)) as Box<dyn SmoothingStrategy>)
            }
            "LGT" => Some(Box::new(SGTSmoother::new()) as Box<dyn SmoothingStrategy>),
            _ => None,
        }
    };

    let smoothing_strategy = smoothing_strategy.unwrap();

    info!(log, "Loading quants from {:?}", &app_info.quant_path);
    let quant_map = load_express_quants(&app_info.quant_path).unwrap();
    let txs_w_mapped_reads = load_express_txps_w_mapped_reads(&app_info.quant_path).unwrap();

    info!(
        log,
        "Loading held-out ECs from {:?}", &app_info.val_ecs_path
    );
    let te_ecs = EqClassCollection::from_path(&app_info.val_ecs_path).unwrap();

    let ec_eval_results = perplexity(
        te_ecs,
        &quant_map,
        &txs_w_mapped_reads,
        smoothing_strategy.as_ref(),
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
    txs_w_mapped_reads: &HashSet<String>,
    smoothing_strategy: &dyn SmoothingStrategy,
    _log: &slog::Logger,
) -> EqClassCollectionEval {
    let train_txs = txs_w_mapped_reads;
    let smoothed_quant_map = smoothing_strategy.smooth_quants(quant_map);

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
            let eta_smooth = smoothed_quant_map[t].tpm / 1e6;
            ec_perp_smooth += w * eta_smooth;
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
