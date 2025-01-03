use std::{ops::Mul, vec};

use curve25519_dalek::{traits::Identity, RistrettoPoint, Scalar};

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

pub fn inv_vector(vector: &Vec<Scalar>) -> Vec<Scalar> {
    vector.iter().map(|x| x.invert()).collect()
}

pub fn inner_product(vector_1: &Vec<Scalar>, vector_2: &Vec<RistrettoPoint>) -> RistrettoPoint {
    let mut value = vector_2[0];

    for (a, b) in vector_1.iter().zip(vector_2.iter()) {
        value = value + (a * b);
    }
    value = value - vector_2[0];

    value
}

pub fn diagonal_ss_sum(vector1: &mut Vec<Scalar>, vector2: &mut Vec<Scalar> ) -> Scalar {
    if vector1.len() % 2 !=0 {
        vector1.insert(0, Scalar::from(0u8));
        vector2.insert(0, Scalar::from(0u8));
    }
    let vec1:Vec<Scalar> = vector1.iter().enumerate().filter(|(i,_)| i % 2 == 0).map(|(_, &x)| x).collect();
    let vec2:Vec<Scalar> = vector2.iter().enumerate().filter(|(i,_)| i % 2 != 0).map(|(_, &x)| x).collect();
    vec1.iter().zip(vec2.iter()).map(|(x, y)| x * y).sum()
}

pub fn diagonal_vs_sum(vector1: &mut Vec<RistrettoPoint>, vector2: &mut Vec<Scalar>) -> RistrettoPoint {
    if vector1.len() % 2 !=0 {
        vector1.insert(0, RistrettoPoint::identity());
        vector2.insert(0, Scalar::from(0u8));
    }
    let vec1:Vec<RistrettoPoint> = vector1.iter().enumerate().filter(|(i,_)| i % 2 == 0).map(|(_, &x)| x).collect();
    let vec2:Vec<Scalar> = vector2.iter().enumerate().filter(|(i,_)| i % 2 != 0).map(|(_, &x)| x).collect();
    vec1.iter().zip(vec2.iter()).map(|(x, y)| x * y).sum()
}

pub fn diagonal_sv_sum(vector1: &mut Vec<Scalar>, vector2: &mut Vec<RistrettoPoint> ) -> RistrettoPoint {
    if vector1.len() % 2 !=0 {
        vector1.insert(0, Scalar::from(0u8));
        vector2.insert(0, RistrettoPoint::identity());
    }
    let vec1:Vec<Scalar> = vector1.iter().enumerate().filter(|(i,_)| i % 2 == 0).map(|(_, &x)| x).collect();
    let vec2:Vec<RistrettoPoint> = vector2.iter().enumerate().filter(|(i,_)| i % 2 != 0).map(|(_, &x)| x).collect();
    vec1.iter().zip(vec2.iter()).map(|(x, y)| x * y).sum()
}

pub fn vector_sub(vector1: &Vec<Scalar>, vector2: &Vec<Scalar>) -> Vec<Scalar> {
    vector1.iter().zip(vector2.into_iter()).map(|(x, y)| x - y).collect()
}

pub fn vector_add(vector1: &Vec<Scalar>, vector2: &Vec<Scalar>) -> Vec<Scalar> {
    vector1.iter().zip(vector2.into_iter()).map(|(x, y)| x + y).collect()
}

pub fn hadamard_multiply(vector1: &Vec<Scalar>, vector2: &Vec<Scalar>) -> Vec<Scalar> {
    vector1.iter().zip(vector2.into_iter()).map(|(x, y)| x * y).collect()
}

pub fn points_hadamard_multiply(vector1: &Vec<Scalar>, vector2: &Vec<RistrettoPoint>) -> Vec<RistrettoPoint> {
    vector1.iter().zip(vector2.into_iter()).map(|(x, y)| x * y).collect()
}

pub fn vec_scalar_mul(vector: &Vec<Scalar>, scalar: &Scalar) -> Vec<Scalar> {
    vector.iter().map(|x| x * scalar).collect()
}