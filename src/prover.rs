use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand::rngs::OsRng;
use crate::operations::to_bin;
use crate::generator::{gen_basis_vectors, gen_v_points};

pub struct ASVcommitment {
    commit_a: RistrettoPoint,
    commit_s: RistrettoPoint,
    commit_v: RistrettoPoint,
}

impl ASVcommitment {
    pub fn asv_commits(v: u64, range: usize) -> Self {
        let seed = b"G and H basis seed";
        let g_basis = gen_basis_vectors(range, seed, "g_basis");
        let h_basis = gen_basis_vectors(range, seed, "h_basis");

        let rng = OsRng;
        let points = gen_v_points(rng);
        let g_i = points[0];
        let h_i = points[1];

        let [al, ar] = ASVcommitment::al_ar(v, range);
    }

    pub fn al_ar(v: u64, range: usize) -> [Vec<i64>; 2] {
        let mut al= to_bin(v.clone());
        while al.len() < range {
            al.push(0)
        }

        let ar = al.iter().map(|&x| x - 1).collect();

        [al, ar]
    }
}