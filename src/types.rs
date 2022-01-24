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
    pub vals: Vec<i32>,
}

impl WaveFront {
    pub fn new(hi: i32, lo: i32) -> Self {
        let len = num::abs_sub(hi, lo) as usize + 1;

        Self {
            hi,
            lo,
            vals: vec![-10; len],
        }
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
    pub wavefront_set: Vec<WaveFrontSet>,

		pub min_k: i32,
		pub max_k: i32,
		pub a_k: usize
}

impl WaveFronts {
		/// The scores should always be positive numbers
		pub fn get(&self, score: usize) -> &WaveFrontSet {
				&self.wavefront_set[score]
		}

		pub fn option_get(&self, score: usize) -> Option<&WaveFrontSet> {
				self.wavefront_set.get(score)
		}


		// Allocate wavefronts with given score
		pub fn allocate_wavefronts(
				&mut self,
				score: u32,
				lo: i32,
				hi: i32
		) -> Result<(), &str>
		{
				// should only add what is necessary
				let prev_score = self.wavefront_set.len() - 1;

				if prev_score >= score as usize {
						// eprintln!("previous score {} score {}", prev_score, score);
						return Err("[types::allocate_wavefronts] fuckery detected");
				}

				for _ in prev_score..(score as usize) {
						let wf_set = WaveFrontSet {
								i: WaveFront::new(hi, lo),
								d: WaveFront::new(hi, lo),
								m: WaveFront::new(hi, lo),
						};

						self.wavefront_set.push(wf_set);
				}

				Ok(())
		}
}
