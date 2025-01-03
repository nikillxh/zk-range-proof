use curve25519_dalek::traits::Identity;
use curve25519_dalek::{RistrettoPoint, Scalar};

use crate::generator::GlobalPoints;
use crate::prover::BulletProof;
use crate::verifier::{Generatives, BulletVerify};

pub fn prove_commitments_log(count: usize, points: &GlobalPoints, prover: &mut BulletProof, verifier: &mut BulletVerify, gen: &Generatives) -> () {
    if prover.a().len() <= 1 {
        verifier.verify([prover.a(), prover.b()], count, points, gen);
    } else {
        let [left, right] = BulletProof::compute_diagonal([prover.a().clone(), prover.b().clone()], [&mut prover.g_basis(), &mut prover.h_basis()], points.g_i());
        verifier.u_gen();
        let u_random = verifier.u_random();
        verifier.compute([left, right]);
        prover.update_diagonals([left, right]);
        prover.compute(u_random);
        prove_commitments_log(count, points, prover, verifier, gen);
    }
}

pub fn fold_scalar(a: &mut Vec<Scalar>, u: Scalar) -> Vec<Scalar> {
    if a.len() % 2 !=0 {
        a.insert(0, Scalar::from(0u8));
    }
    let folded: Vec<Scalar> = a.chunks(2).map(|chunk| {
        match chunk {
            [x, y] => (x * u) + (y * u.invert()),
            _ => panic!("Scalar folding error")
        }
    }).collect();

    // assert_eq!(a.len(), folded.len() * 2, "Scalar folded incorrectly");
    folded
}

pub fn fold_vector(a: &mut Vec<RistrettoPoint>, u: Scalar) -> Vec<RistrettoPoint> {
    if a.len() % 2 !=0 {
        a.insert(0, RistrettoPoint::identity());
    }

    let folded: Vec<RistrettoPoint> = a.chunks(2).map(|chunk| {
        match chunk {
            [x, y] => (x * u) + (y * u.invert()),
            _ => panic!("Scalar folding error")
        }
    }).collect();

    // assert_eq!(a.len(), folded.len() * 2, "Vector folded incorrectly");
    folded
}

