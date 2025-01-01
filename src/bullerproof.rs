use curve25519_dalek::{RistrettoPoint, Scalar};
use crate::{generator::GlobalPoints, prover::BulletProof, verifier::{self, BulletVerify}};

pub fn prove_commitments_log(commit_c: RistrettoPoint, points: GlobalPoints, prover: BulletProof, verifier: BulletVerify) -> bool {
    if prover.a().len() == 1 {

    } else
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



