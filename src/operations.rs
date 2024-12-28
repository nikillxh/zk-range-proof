use std::vec;

use curve25519_dalek::{RistrettoPoint, Scalar};

pub fn to_bin(mut v: u64) -> Vec<i64> {
    let mut al: Vec<i64> = vec![];

    while v > 0 {
        al.push((v % 2) as i64);
        v /= 2;
    }

    al
}

pub fn to_dec(al: &Vec<i64>) -> u64 {
    let v: i64 = al.iter().enumerate().map(|(i, &x)| 2i64.pow(i as u32) * x).sum();

    v as u64
}

pub fn scalarize(vector: &mut Vec<i64>) -> Vec<Scalar> {
    fn handle_val(value: i64) -> Scalar {
        if value >= 0 {
            Scalar::from(value as u64)
        } else {
            Scalar::from(0u64) - Scalar::from(value.abs() as u64)
        }
    }

    let scalar_vector = vector.iter_mut().map(|&mut x| handle_val(x)).collect();

    scalar_vector
}

pub fn inner_product(vector_1: &Vec<Scalar>, vector_2: &Vec<RistrettoPoint>) -> RistrettoPoint {
    let mut value = vector_2[0];

    for (a, b) in vector_1.iter().zip(vector_2.iter()) {
        value = value + (a * b);
    }
    value = value - vector_2[0];

    value
}