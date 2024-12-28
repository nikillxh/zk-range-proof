use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand::rngs::OsRng;
use crate::operations::{scalarize, inner_product, to_bin};
use crate::generator::{gen_basis_vectors, gen_scalars, gen_v_points};

pub struct ASVcommitment {
    commit_a: RistrettoPoint,
    commit_s: RistrettoPoint,
    commit_v: RistrettoPoint,
}

impl ASVcommitment {
    pub fn asv_commits(v: u64, range: usize) -> Self {
        let seed = b"G and H basis seed";
        let mut g_basis = gen_basis_vectors(range, seed, "g_basis");
        let mut h_basis = gen_basis_vectors(range, seed, "h_basis");

        let mut rng = OsRng;
        let points = gen_v_points(rng);
        let g_i = points[0];
        let b_i = points[1];

        let salt_alpha = Scalar::random(&mut rng);
        let salt_beta = Scalar::random(&mut rng);
        let salt_gamma = Scalar::random(&mut rng);

        let [mut al, mut ar] = ASVcommitment::al_ar(v, range);
        let sl = gen_scalars(range);
        let sr = gen_scalars(range);

        println!("ASV commitment pre-requirements completed.");

        let commit_a_val = inner_product(&scalarize(&mut al), &mut g_basis)
            + inner_product(&scalarize(&mut ar), &mut h_basis) + (salt_alpha * b_i);

        let commit_s_val = inner_product(&sl, &mut g_basis)
            + inner_product(&sr, &mut h_basis) + (salt_beta * b_i);
        
        let commit_v_val = (Scalar::from(v) * g_i) + (salt_gamma * b_i);

        println!("ASV commitments sent to verifier!");

        Self {
            commit_a: commit_a_val,
            commit_s: commit_s_val,
            commit_v: commit_v_val,
        }
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