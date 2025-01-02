use curve25519_dalek::{RistrettoPoint, Scalar};

use crate::generator::GlobalPoints;
use crate::prover::BulletProof;
use crate::verifier::{Generatives, BulletVerify};

pub fn prove_commitments_log(count: usize, points: GlobalPoints, mut prover: BulletProof, mut verifier: BulletVerify, gen: Generatives) -> () {
    if prover.a().len() <= 1 {
        verifier.verify([prover.a(), prover.b()], count, points, gen);
    } else {
        let [left, right] = BulletProof::compute_diagonal([prover.a().clone(), prover.b().clone()], &points);
        let u_random = verifier.u_update();
        verifier.compute([left, right]);
        prover.compute(u_random);
        prover.a_fold();
        prover.b_fold();
        prove_commitments_log(count, points, prover, verifier, gen);
    }
}

pub fn fold_scalar(mut a: Vec<Scalar>, u: Scalar) -> Vec<Scalar> {
    if a.len() % 2 !=0 {
        a.insert(0, Scalar::from(0u8));
    }
    let folded = a.chunks(2).map(|chunk| {
        match chunk {
            [x, y] => (x * u) + (y * u.invert()),
            _ => panic!("Scalar folding error")
        }
    }).collect();

    folded
}

pub fn fold_vector(mut a: Vec<RistrettoPoint>, u: Scalar) -> Vec<RistrettoPoint> {
    let mut folded: Vec<RistrettoPoint> = vec![];

    if a.len() % 2 !=0 {
        folded.push(a.remove(0) * u.invert());
    }

    folded.extend(a.chunks(2).map(|chunk| {
        match chunk {
            [x, y] => (x * u) + (y * u.invert()),
            _ => panic!("Vector folding error")
        }
    }));

    folded
}



