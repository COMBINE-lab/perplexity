// Smoothing strategies for perplexity computation
use seine::salmon::QuantEntry;
use std::collections::{BTreeMap, HashMap};
use std::default::Default;

pub trait SmoothingStrategy {
    fn smooth_quants(&self, quant_map: &HashMap<String, QuantEntry>)
        -> HashMap<String, QuantEntry>;
}

pub struct PerTPMSmoother {
    pub alpha: f64,
}

pub struct PerTxpSmoother {
    pub alpha: f64,
}

pub struct SGTSmoother {}

impl Default for PerTxpSmoother {
    fn default() -> Self {
        Self { alpha: 1.0 }
    }
}

impl Default for PerTPMSmoother {
    fn default() -> Self {
        Self { alpha: 1e-8 }
    }
}

impl PerTxpSmoother {
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl PerTPMSmoother {
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl SGTSmoother {
    pub fn new() -> Self {
        Self {}
    }
}

impl SmoothingStrategy for PerTxpSmoother {
    fn smooth_quants(
        &self,
        quant_map: &HashMap<String, QuantEntry>,
    ) -> HashMap<String, QuantEntry> {
        let mut smoothed_quants = HashMap::new();

        // compute the new denominator
        let mut denom = 0.0;
        for (_, v) in quant_map.iter() {
            denom += (v.num_reads + self.alpha) / v.efflen;
        }

        for (k, v) in quant_map.iter() {
            // add $\alpha$ reads per transcript
            let num_reads = v.num_reads + self.alpha;
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
}

impl SmoothingStrategy for PerTPMSmoother {
    fn smooth_quants(
        &self,
        quant_map: &HashMap<String, QuantEntry>,
    ) -> HashMap<String, QuantEntry> {
        let mut smoothed_quants = HashMap::new();

        let m = quant_map.len();

        let denom = 1. + (m as f64 * self.alpha);

        // reconstituted counts
        let mut total_reads = 0.;
        for (_, v) in quant_map.iter() {
            total_reads = total_reads + v.num_reads;
        }

        for (k, v) in quant_map.iter() {
            let eta = ((v.tpm / 1e6) + self.alpha) / denom;
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
}

impl SmoothingStrategy for SGTSmoother {
    fn smooth_quants(
        &self,
        quant_map: &HashMap<String, QuantEntry>,
    ) -> HashMap<String, QuantEntry> {
        let mut smoothed_quants = HashMap::new();

        // 1. Calculate frequencies of frequencies
        let (rs, nrs) = freq_of_freqs(quant_map);
        let zs = calc_zs(&rs, &nrs);

        // 2. Fit log-linear function to denoise freqs of freqs
        let (slope, b) = {
            let rs = rs.iter().map(|r| *r as f64).collect();
            log_linfit(&rs, &zs)
        };

        // 3. Iterate over input frequencies to determine when to
        //    switch from turing-estimates to LGT-estimates

        let mut use_lgt = false;
        let m = rs.len();
        let mut rstar = vec![-1.; m];
        for (i, r) in rs.iter().enumerate() {
            let nr_plus1 = if (i < (m - 1)) && ((r + 1) == rs[i + 1]) {
                nrs[i + 1] as f64
            } else {
                0.
            };
            let nr = nrs[i] as f64;
            let r = *r as f64;

            let sr_plus1 = (b + slope * (r + 1.).ln()).exp();
            let sr = (b + slope * r.ln()).exp();

            let lgt_r = (r + 1.) * sr_plus1 / sr;
            let t_r = (r + 1.) * nr_plus1 / nr;

            if !use_lgt {
                // check variance.
                let mut var = (r + 1.).powf(2.);
                var *= nr_plus1 / nr.powf(2.);
                var *= 1. + nr_plus1 / nr;
                use_lgt = (lgt_r - t_r).abs() < (var * 1.65)
            }

            if use_lgt {
                rstar[i] = lgt_r;
            } else {
                rstar[i] = t_r;
            }
        }

        // 4. Compute total probability mass for unseen txps
        let n_tot: f64 = rs
            .iter()
            .zip(nrs.iter())
            .map(|(r, nr)| (*r as f64) * (*nr as f64))
            .sum();

        let p0 = nrs[0] as f64 / n_tot;

        // 5. Compute adjusted total observations
        let big_n_prime: f64 = rstar
            .iter()
            .zip(nrs.iter())
            .map(|(r, nr)| r * (*nr as f64))
            .sum();

        let n_hat = big_n_prime / (1. - p0);

        let n_0 = n_hat - big_n_prime;

        // 6. uniformly reweight p0 proportional to efflens of unseen txps
        let total_unexpressed_nucs = {
            let mut sum = 0.;
            for (_, v) in quant_map.iter() {
                if v.num_reads == 0. {
                    sum += v.efflen;
                }
            }
            sum
        };

        let per_nuc = n_0 / total_unexpressed_nucs;

        let mut rstar_map = HashMap::new();
        for (r, r_adj) in rs.iter().zip(rstar.iter()) {
            rstar_map.insert(*r, *r_adj);
        }

        // 7. Compute new denominator for TPM calculation
        let mut denom = 0.;
        for (_, v) in quant_map.iter() {
            let r = v.num_reads.round() as usize;
            let num_reads = if r == 0 {
                v.efflen * per_nuc
            } else {
                *rstar_map.get(&r).expect(&format!("{:?}", v))
            };
            denom += num_reads / v.efflen;
        }

        let denom = denom;

        // 8. Adjust TPMs and expected no. reads
        for (k, v) in quant_map.iter() {
            let r = v.num_reads.round() as usize;
            let num_reads = if r == 0 {
                v.efflen * per_nuc
            } else {
                *rstar_map.get(&r).expect(&format!("{:?}", v))
            };

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
}

// Good turing smoothing helpers
fn freq_of_freqs(quant_map: &HashMap<String, QuantEntry>) -> (Vec<usize>, Vec<usize>) {
    let mut freqs = BTreeMap::new();
    for (_, v) in quant_map.iter() {
        let num_reads = v.num_reads;
        let num_reads = num_reads.round() as usize;

        let freq = freqs.entry(num_reads).or_insert(0);
        *freq += 1;
    }

    let mut rs: Vec<usize> = freqs.keys().cloned().collect();
    let mut nrs: Vec<usize> = freqs.values().cloned().collect();

    // silly sanity checks...
    assert!(rs[0] == 0);
    assert!(rs[1] == 1);
    assert!(nrs[1] > 0);

    rs.remove(0);
    nrs.remove(0);

    (rs, nrs)
}

fn mean_f64(numbers: &Vec<f64>) -> f64 {
    let sum: f64 = numbers.iter().sum();
    sum / (numbers.len() as f64)
}

fn calc_zs(rs: &Vec<usize>, nrs: &Vec<usize>) -> Vec<f64> {
    let m = rs.len();
    let mut zs = vec![0.; m];

    for (i, r) in rs.iter().enumerate() {
        let q = if i == 0 { 0 } else { rs[i - 1] };
        let t = if i == (m - 1) { 2 * r - q } else { rs[i + 1] };

        let nr = nrs[i] as f64;
        let denom = (t - q) as f64;

        let z = 2. * nr / denom;

        zs[i] = z;
    }
    zs
}

fn log_linfit(x: &Vec<f64>, y: &Vec<f64>) -> (f64, f64) {
    // Use (X.T X)^{-1}   X.T Y
    // zero mean X and y

    let x: Vec<f64> = x.iter().map(|x| x.ln()).collect();
    let y: Vec<f64> = y.iter().map(|y| y.ln()).collect();

    let mean_x = mean_f64(&x);
    let mean_y = mean_f64(&y);

    let x: Vec<f64> = x.iter().map(|x| x - mean_x).collect();
    let y: Vec<f64> = y.iter().map(|y| y - mean_y).collect();

    let xy = x.iter().zip(y.iter());

    let top: f64 = xy.into_iter().map(|(x, y)| x * y).sum();
    let bottom: f64 = x.iter().map(|x| x * x).sum();

    let slope = top / bottom;
    let bias = mean_y - slope * mean_x;
    (slope, bias)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init() {
        let smooth = PerTPMSmoother::new(1.);
        assert_eq!(smooth.alpha, 1.);

        let smooth = PerTxpSmoother::new(5.);
        assert_eq!(smooth.alpha, 5.);

        let smooth = PerTPMSmoother::default();
        assert_eq!(smooth.alpha, 1e-8);

        let smooth = PerTxpSmoother::default();
        assert_eq!(smooth.alpha, 1.);
    }
}
