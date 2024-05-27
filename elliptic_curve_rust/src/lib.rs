use num_bigint::BigUint;
use std::iter::Product;

#[derive(PartialEq, Clone, Debug)]
enum Point {
    Coor(BigUint, BigUint),
    Identity,
}
struct EllipticCurve {
    // y^2 = x^2 + a * x + b
    a: BigUint,
    b: BigUint,
    p: BigUint,
}

impl EllipticCurve {
    fn add(&self, c: &Point, d: &Point) -> Point {
        assert!(self.is_on_curve(c), "Point is not on the curve");
        assert!(self.is_on_curve(d), "Point is not on the curve");
        assert!(*c != *d, "Points should not be the same");

        match (c, d) {
            (Point::Identity, _) => d.clone(),
            (_, Point::Identity) => c.clone(),
            (Point::Coor(x1, y1), Point::Coor(x2, y2)) => {
                let y1plusy2 = FiniteField::add(&y1, &y2, &self.p);
                if x1 == x2 && y1plusy2 == BigUint::from(0u32){
                    return Point::Identity;
                }
                // s = (y2 - y1) / (x2 - x1) mod p
                // x3 = s^2 - x1 - x2 mod p
                // y3 = s(x1 - x3) - y1 mod p
                let y2minusy1 = FiniteField::subtract(y2, y1, &self.p);
                let x2minusx1 = FiniteField::subtract(x2, x1, &self.p);
                let s = FiniteField::divide(&y2minusy1, &x2minusx1, &self.p);
                let s2 = s.modpow(&BigUint::from(2u32), &self.p);
                let s2minusx1 = FiniteField::subtract(&s2, x1, &self.p);
                let x3 = FiniteField::subtract(&s2minusx1, x2, &self.p);
                let x1minusx3 = FiniteField::subtract(x1, &x3, &self.p);
                let sx1minusx3 = FiniteField::mult(&s, &x1minusx3, &self.p);
                let y3 = FiniteField::subtract(&sx1minusx3, &y1, &self.p);
                Point::Coor(x3, y3)
            }
        }
    }

    fn double(&self, c: &Point) -> Point {
        assert!(self.is_on_curve(c), "Point is not on the curve");

        match c {
            Point::Identity => Point::Identity,
            Point::Coor(x1, y1) => {
                // s = (3* x1^2 + a) / (2 * y1) mod p
                // x3 = s^2 - 2 * x1 mod p
                // y3 = s(x1 - x3) - y1 mod p
                let numerator = x1.modpow(&BigUint::from(2u32), &self.p);
                let numerator = FiniteField::mult(&BigUint::from(3u32), &numerator, &self.p);
                let numerator = FiniteField::add(&self.a, &numerator, &self.p);

                let denominator = FiniteField::mult(&BigUint::from(2u32), y1, &self.p);

                let s = FiniteField::divide(&numerator, &denominator, &self.p);

                let s2 = s.modpow(&BigUint::from(2u32), &self.p);
                let s2minusx1 = FiniteField::subtract(&s2, x1, &self.p);
                let x3 = FiniteField::subtract(&s2minusx1, x1, &self.p);
                let x1minusx3 = FiniteField::subtract(x1, &x3, &self.p);
                let sx1minusx3 = FiniteField::mult(&s, &x1minusx3, &self.p);
                let y3 = FiniteField::subtract(&sx1minusx3, &y1, &self.p);
                Point::Coor(x3, y3)

                
            }
        }
    }

    fn scalar_mul(&self, c: &Point, d: &BigUint) -> Point {
        // addition/doubling algorithm
        // B = d * A
        let mut t = c.clone();
        for i in (0..(d.bits() - 1)).rev() {
            t = self.double(&t);
            if d.bit(i) {
                t = self.add(&t, c);
            }
        }
        t
    }

    fn is_on_curve(&self, c: &Point) -> bool{
        
        match c {
            Point::Coor(x, y) => {
                // y^2 = x^3 + a * x + b
                let y2 = y.modpow(&BigUint::from(2u32), &self.p);
                let x3 = x.modpow(&BigUint::from(3u32), &self.p);
                let ax = FiniteField::mult(&self.a, x, &self.p);
                let x3plusax = FiniteField::add(&x3, &ax, &self.p);
                y2 == FiniteField::add(&x3plusax, &self.b, &self.p)
            },
            Point::Identity => true,
        }
    }

}


struct FiniteField {}

impl FiniteField{
    fn add(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint{
        // c + d = r mod p
        let r = c + d;
        r.modpow(&BigUint::from(1u32), p)
    }

    fn mult(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint{
        // c * d = r mod p
        let r = c * d;
        r.modpow(&BigUint::from(1u32), p)
    }

    fn inv_addition(c: &BigUint, p: &BigUint) -> BigUint{
        // -c mod p
        assert!(c < p, "number is >= p");
        p - c
    }

    fn subtract(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint{
        let d_inv = FiniteField::inv_addition(d, p);
        FiniteField::add(c, &d_inv, p)
    }

    fn inv_multiplication(c: &BigUint, p: &BigUint) -> BigUint{
        // only for p prime
        // c^(-1) mod p = c^(p-2) mod p
        c.modpow(&(p - BigUint::from(2u32)), p)
    }

    fn divide(c: &BigUint, d: &BigUint, p: &BigUint) -> BigUint{
        let d_inv = FiniteField::inv_multiplication(d, p);
        FiniteField::mult(c, &d_inv, p)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_add(){
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::add(&c, &d, &p);

        assert_eq!(r, BigUint::from(3u32));
    }

    #[test]
    fn test_mul(){
        let c = BigUint::from(4u32);
        let d = BigUint::from(10u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::mult(&c, &d, &p);

        assert_eq!(r, BigUint::from(7u32));
    }

    #[test]
    fn test_inv_addition(){
        let c = BigUint::from(4u32);
        let p = BigUint::from(51u32);

        let r = FiniteField::inv_addition(&c, &p);

        assert_eq!(r, BigUint::from(47u32));
    }

    #[test]
    fn test_inv_multiplication(){
        let c = BigUint::from(4u32);
        let p = BigUint::from(11u32);

        let r = FiniteField::inv_multiplication(&c, &p);

        assert_eq!(r, BigUint::from(3u32));
    }

    #[test]
    fn test_ec_point_addition () {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve{
            a: BigUint::from(2u32), 
            b: BigUint::from(2u32), 
            p: BigUint::from(17u32)
        };

        // (6,3) + (5,1) = (10, 6)
        let p1 = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));
        let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(10u32), BigUint::from(6u32));

        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);

    }

    #[test]
    fn test_ec_point_addition_opposite () {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve{
            a: BigUint::from(2u32), 
            b: BigUint::from(2u32), 
            p: BigUint::from(17u32)
        };

        // (5,16) + (5,1) = Point::Identity
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(16u32));
        let p2 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Identity;

        let res = ec.add(&p1, &p2);
        assert_eq!(res, pr);

    }

    #[test]
    fn test_ec_point_doubeling () {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve{
            a: BigUint::from(2u32), 
            b: BigUint::from(2u32), 
            p: BigUint::from(17u32)
        };

        // (5,1) + (5,1) = (6,3)
        let p1 = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));
        let pr = Point::Coor(BigUint::from(6u32), BigUint::from(3u32));

        let res = ec.double(&p1);
        assert_eq!(res, pr);

    }

    #[test]
    fn test_ec_scalar_multiplication () {
        // y^2 = x^3 + 2x + 2 mod 17  |G| = 19 19 * A = I
        let ec = EllipticCurve{
            a: BigUint::from(2u32), 
            b: BigUint::from(2u32), 
            p: BigUint::from(17u32)
        };

        
        let c = Point::Coor(BigUint::from(5u32), BigUint::from(1u32));

        // 10 * (5,1) = (7,11)

        let pr = Point::Coor(BigUint::from(7u32), BigUint::from(11u32));

        let res = ec.scalar_mul(&c, &BigUint::from(10u32));
        assert_eq!(res, pr);

        // 19 * (5,1) = I

        let pr = Point::Identity;

        let res = ec.scalar_mul(&c, &BigUint::from(19u32));
        assert_eq!(res, pr);

    }
}

