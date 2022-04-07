/*!
Configs and related WFA types
 */

use crate::utils;

// ----------------------
//         Config
// ----------------------
// TODO: use u8
pub struct Penalties {
    pub mismatch: i32,
    pub matches: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

pub struct Config {
    pub adapt: bool,
    // pub segment_length: u32, // segment size in bytes
    // pub step_size: u32,
    // pub thread_count: usize,
    pub verbosity: u8,
    pub penalties: Penalties,
}

// ----------------------
//     Core types
// ----------------------
#[derive(Debug, PartialEq, Eq)]
pub enum WfType {
    D,
    I,
    M,
}

#[derive(Debug, PartialEq, Eq)]
pub enum BacktraceOperation {
    MatchMismatch,
    Insertion,
    Deletion,
}

/// The a single wavefront with a score
/// The furthest reaching point of a single wavefront
#[derive(Debug, Clone)]
pub struct WaveFront {
    /// the highest diagonal touched by the wavefront
    pub hi: i32,

    /// the lowest diagonal touched by the wavefront
    pub lo: i32,

    /// The offsets of each diagonal between hi and lo
    /// vals\[0\] is the score of the wavefront at diagonal hi and
    /// vals\[<last>\] is the score of the wavefront at diagonal lo
    /// length is (hi - lo) + 1
    pub offsets: Vec<i32>,
}

impl WaveFront {
    pub fn new(hi: i32, lo: i32) -> Self {
        let len = utils::new_compute_wave_length(lo, hi);

        Self {
            hi,
            lo,
            offsets: vec![0; len],
        }
    }

    pub fn len(&self) -> usize {
        // TODO merge with utils

        let wave_len = utils::new_compute_wave_length(self.lo, self.hi);
        // TODO: remove
        assert_eq!(wave_len, self.offsets.len());

        wave_len
    }

    pub fn k_index(&self, k: i32) -> usize {
        utils::new_compute_k_index(k, self.lo, self.hi)
    }

    pub fn in_bounds(&self, k: i32) -> bool {
        utils::k_in_bounds(k, self.lo, self.hi)
    }
    /// return the offset at the k diagonal
    /// Computes the k-index internally
    /// takes the diagonal k (not k-index)
    pub fn get_offset(&self, k: i32) -> Option<&i32> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);
        self.offsets.get(k_index)
    }

    pub fn get_offset_mut(&mut self, k: i32) -> Option<&mut i32> {
        if !utils::k_in_bounds(k, self.lo, self.hi) {
            // eprintln!("out of bounds mut");
            return None;
        }

        let k_index = utils::new_compute_k_index(k, self.lo, self.hi);
        self.offsets.get_mut(k_index)
    }
}

/// The set of wavefronts at a certain score
#[derive(Debug, Clone)]
pub struct WaveFrontSet {
    /// insertion wavefront
    pub i: Option<WaveFront>,

    /// deletion wavefront
    pub d: Option<WaveFront>,

    /// match wavefront
    /// $\tilde{M}_{s,k}$ is the value of the m wavefront at diagonal k
    pub m: Option<WaveFront>,
}

/// All the wavefronts
#[derive(Clone)]
pub struct WaveFronts {
    /// The set of wavefronts with each score, the index represents the score
    /// and, each element is a wavefront.
    /// WF_s is wavefront_set\[s\]
    pub wavefront_set: Vec<Option<WaveFrontSet>>,

    pub min_k: isize, // -qlen
    pub max_k: isize, // tlen
    pub a_k: i32,
}

impl WaveFronts {
    /// The scores should always be positive numbers
    pub fn get(&self, score: usize) -> &Option<WaveFrontSet> {
        &self.wavefront_set[score]
    }

    pub fn option_get(&self, score: usize) -> Option<&WaveFrontSet> {
        self.wavefront_set
            .get(score)
            .and_then(|maybe_wf_set| maybe_wf_set.as_ref())
    }

    pub fn get_wavefronts(&self, score: usize) -> Option<&WaveFrontSet> {
        self.wavefront_set
            .get(score)
            .and_then(|maybe_wf_set| maybe_wf_set.as_ref())
    }

    pub fn get_m_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.m.as_ref(),
            _ => None,
        }
    }

    pub fn get_i_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.i.as_ref(),
            _ => None,
        }
    }

    pub fn get_d_wavefront(&self, score: i32) -> Option<&WaveFront> {
        if score < 0 {
            return None;
        }

        let score = score as usize;
        let maybe_wf_set: Option<&WaveFrontSet> = self.option_get(score);
        match maybe_wf_set {
            Some(v) => v.d.as_ref(),
            _ => None,
        }
    }

    pub fn set_i_d_m(&mut self, score: usize) {
        let wf_set: &mut Option<WaveFrontSet> = self.wavefront_set.get_mut(score).unwrap();
        let wf_set: &mut Option<WaveFrontSet> = &mut self.wavefront_set[score];
    }

    pub fn len(&self) -> usize {
        self.wavefront_set.len()
    }

    pub fn max_score(&self) -> u32 {
        (self.len() - 1) as u32
    }

    // Allocate wavefronts (I, D, or M) for the given score
    pub fn allocate_wavefronts(
        &mut self,
        score: u32,
        lo: i32,
        hi: i32,
        wavefronts_to_allocate: &Vec<WfType>,
    ) -> Result<(), &str> {
        // should only add what is necessary
        let max_score = self.max_score();
        let len = num::abs_sub(hi, lo) as usize + 1;

        if max_score >= score {
            // we are trying to add a score that exists
            // eprintln!("previous score {} score {}", prev_score, score);
            return Err("[types::allocate_wavefronts] fuckery detected");
        }

        for index in max_score + 1..=score {
            if index == score {
                let wf_set = WaveFrontSet {
                    i: {
                        if wavefronts_to_allocate.contains(&WfType::I) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                    d: {
                        if wavefronts_to_allocate.contains(&WfType::D) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                    m: {
                        if wavefronts_to_allocate.contains(&WfType::M) {
                            Some(WaveFront::new(hi, lo))
                        } else {
                            None
                        }
                    },
                };

                self.wavefront_set.push(Some(wf_set));
            } else {
                self.wavefront_set.push(None);
            }
        }
        Ok(())
    }
}
