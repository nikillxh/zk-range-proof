use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use operations::{to_bin, vector_sub, scalarize};
use generator::GlobalPoints;
use prover::{ASVcommitment, Polycommitment, Salts, T1T2commitment};
use rand::rngs::OsRng;
use verifier::Generatives;

mod prover;
mod generator;
mod operations;
mod verifier;
mod bullerproof;

fn main() {
    // Prover
    let prover_value: u64 = 54;
    let al = scalarize(&mut to_bin(prover_value));
    let range = al.len();
    let ar = vector_sub( &al, &vec![Scalar::from(1u32)]);
    let salt = Salts::init();

    // Global
    let points= GlobalPoints::gen_global(range);

    // Prover
    let asv = ASVcommitment::compute(prover_value, range, &salt, &points);

    // Verifier
    let gen = Generatives::init(range, &points);
    let [y, z] = [gen.y(), gen.z()];

    // Prover
    let t1t2 = T1T2commitment::init(&salt, asv, y, z, range, &points);
    let [commit_t1, commit_t2] = [t1t2.commit_t1(), t1t2.commit_t2()];

    

}
