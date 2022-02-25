#![cfg(test)]
#![allow(dead_code, unused_imports)] // Quiet down needless warnings

use std::fs;
use std::path;
use crate::types;
use crate::utils::*;

// Re-exports
pub use crate::wflambda::wf_align as wflambda_align;
pub use crate::wfa::wf_align as wfa_align;
pub use pretty_assertions::{assert_eq};

// We need #[cfg(test)] guards above the exports and test dependencies

use ndarray::{Array, Array2};
use ndarray_to_img;

pub static TEST_CONFIG: crate::types::Config = crate::types::Config {
    adapt: false,
    verbosity: 0,
    penalties: crate::types::Penalties {
        mismatch: 4,
        matches: 0,
        gap_open: 6,
        gap_extend: 2,
    },
};

pub fn visualize_all(
    all_wavefronts: &types::WaveFronts,
    a_offset: u32,
    each_wf: &Vec<types::WfType>,
    match_positions: &Vec<(i32, i32, usize)>,
    config: &types::Config,
    score: usize,
) {
    if config.verbosity > 1 {
        eprintln!("[utils::visualize_all]");
    }

    if config.verbosity > 2 {
        eprintln!("\tPopulating matrix");
    }

    let dim = (a_offset as usize + 1, a_offset as usize + 1);
    let x = ndarray::Dim(dim);
    let mut matrix: Array2<Option<i32>> = Array::from_elem(x, None);

    // eprintln!("\t\tk\tscore\toffset\t(v,h)");
    for s in (0..=score).rev() {
        let wf_specific: &types::WaveFront = match all_wavefronts.get_m_wavefront(s as i32) {
            Some(m) => m,
            _ => continue,
        };

        let lo = wf_specific.lo;
        let hi = wf_specific.hi;
        let offsets = &wf_specific.offsets;
        let len = offsets.len();

        // eprintln!("\t\tscore: {} lo {} hi {}", s, lo, hi);

        for k in lo..=hi {
            let k_index: usize = compute_k_index(len, k, hi);
            let m_s_k: i32 = offsets[k_index];

            for offset in m_s_k..=m_s_k {
                let v = compute_v(offset, k);
                let h = compute_h(offset, k);

                // eprintln!("offset: {}\tk: {}\tscore: {}\t({}, {})", m_s_k,  k, score, v, h);

                // eprintln!("\t\t{}\t{}\t{}\t({},{})", k, s, offset, v, h);

                if v < 0 || h < 0 || v >= dim.0 as i32 || h >= dim.0 as i32 {
                    continue;
                }

                let v = v as usize;
                let h = h as usize;

                if config.verbosity > 5 {
                    // eprintln!("\t\t({},{})\t{}\t{}", v, h, s, offset);
                }

                matrix[[v, h]] = Some(offset);
            }
        }

        // eprintln!();
    }

    for (v, h, score) in match_positions {
        let v = *v;
        let h = *h;

        if v < 0 || h < 0 || v >= dim.0 as i32 || h >= dim.0 as i32 {
            continue;
        }

        let v = v as usize;
        let h = h as usize;

        if v == h {
            matrix[[v, h]] = Some(*score as i32);
        }
    }

    gen_image(&matrix, config);
}

fn gen_image(matrix: &Array2<Option<i32>>, config: &types::Config) {
    if config.verbosity > 1 {
        eprintln!("[utils::gen_image]");
    }

    let config = ndarray_to_img::Config {
        verbosity: 0,
        with_color: true,
        annotate_image: true,
        draw_diagonal: true,
        draw_boundaries: true,
        scaling_factor: 10,
    };

    let scaled_matrix = ndarray_to_img::scale_matrix(&matrix, &config);

    let out_dir = path::Path::new("./debug/out/");
    // assume this will never fail
    if !out_dir.exists() {
        fs::create_dir_all(out_dir).unwrap();
    }

    let image_path_buf: path::PathBuf = out_dir.join("all.png");
    let image_name: &str = image_path_buf.to_str().unwrap();
    ndarray_to_img::generate_image(&scaled_matrix, &config, image_name).unwrap();
}
