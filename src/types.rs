/*!
wavefront types
 */

use num;


#[derive(Debug)]
pub enum WfType {
		D,
		I,
		M
}

#[derive(PartialEq, Eq)]
pub enum BacktraceOperation {
		MatchMismatch,
		Insertion,
		Deletion
}

/// The a single wavefront with a score
/// The furthest reaching point of a single wavefront
pub struct WaveFront {
    /// the highest diagonal touched by the wavefront
    pub hi: i32,

    /// the lowest diagonal touched by the wavefront
    pub lo: i32,

    /// The offsets of each diagonal between hi and lo
    /// vals\[0\] is the score of the wavefront at diagonal hi and
    /// vals\[<last>\] is the score of the wavefront at diagonal lo
    /// length is (hi - lo) + 1
		// TODO: rename vals to offsets
    pub offsets: Vec<i32>,
}

impl WaveFront {
    pub fn new(hi: i32, lo: i32) -> Self {
        let len = num::abs_sub(hi, lo) as usize + 1;

        Self {
            hi,
            lo,
            offsets: vec![0; len],
        }
    }

		pub fn len(&self) -> usize {
				let wave_length = num::abs_sub(self.hi, self.lo) as usize + 1;
				let offset_len = self.offsets.len();
				if wave_length != offset_len {
						panic!("len error")
				}

				offset_len
		}
}

/// The set of wavefronts at a certain score
pub struct WaveFrontSet {
    /// insertion wavefront
    pub i: WaveFront,

		/// deletion wavefront
    pub d: WaveFront,

    /// match wavefront
    /// $\tilde{M}_{s,k}$ is the value of the m wavefront at diagonal k
    pub m: WaveFront,
}

/// All the wavefronts
pub struct WaveFronts {
    /// The set of wavefronts with each score, the index represents the score
    /// and, each element is a wavefront.
    /// WF_s is wavefront_set\[s\]
    pub wavefront_set: Vec<Option<WaveFrontSet>>,

		pub min_k: i32,
		pub max_k: i32,
		pub a_k: usize
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

		pub fn get_m_wavefront(&self, score: usize) -> Option<&WaveFront> {
				self.wavefront_set
						.get(score)
						.and_then(|maybe_wf_set| match maybe_wf_set {
								Some(wf_set) => Some(&wf_set.m),
								_ => None,
						})
		}

		pub fn get_i_wavefront(&self, score: usize) -> Option<&WaveFront> {
				self.wavefront_set
						.get(score)
						.and_then(|maybe_wf_set| match maybe_wf_set {
								Some(wf_set) => Some(&wf_set.i),
								_ => None,
						})
		}

		pub fn get_d_wavefront(&self, score: usize) -> Option<&WaveFront> {
				self.wavefront_set
						.get(score)
						.and_then(|maybe_wf_set| match maybe_wf_set {
								Some(wf_set) => Some(&wf_set.d),
								_ => None,
						})
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

		// Allocate wavefronts (I, D, & M) for the given score
		pub fn allocate_wavefronts(
				&mut self,
				score: u32,
				lo: i32,
				hi: i32
		) -> Result<(), &str>
		{
				// should only add what is necessary
				let max_score = self.max_score();

				if max_score >= score {
						// we are trying to add a score that exists
						// eprintln!("previous score {} score {}", prev_score, score);
						return Err("[types::allocate_wavefronts] fuckery detected");
				}

				for index in max_score+1..=score {
						if index == score {
								let wf_set = WaveFrontSet {
										i: WaveFront::new(hi, lo),
										d: WaveFront::new(hi, lo),
										m: WaveFront::new(hi, lo),
								};

								self.wavefront_set.push(Some(wf_set));
						} else {
								self.wavefront_set.push(None);
						}
				}

				Ok(())
		}
}
