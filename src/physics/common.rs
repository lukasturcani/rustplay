pub fn vector_sub(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 - x2, y1 - y2, z1 - z2)
}

pub fn vector_sub_(
    (x1, y1, z1): (f32, f32, f32),
    (x2, y2, z2): (f32, f32, f32),
) -> (f32, f32, f32) {
    vector_sub(x1, y1, z1, x2, y2, z2)
}

pub fn vector_add(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 + x2, y1 + y2, z1 + z2)
}

pub fn vector_add_(
    (x1, y1, z1): (f32, f32, f32),
    (x2, y2, z2): (f32, f32, f32),
) -> (f32, f32, f32) {
    vector_add(x1, y1, z1, x2, y2, z2)
}

pub fn vector_mul(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> (f32, f32, f32) {
    (x1 * x2, y1 * y2, z1 * z2)
}

pub fn vector_scalar_mul(scalar: f32, x: f32, y: f32, z: f32) -> (f32, f32, f32) {
    (scalar * x, scalar * y, scalar * z)
}

pub fn vector_sum_((x, y, z): (f32, f32, f32)) -> f32 {
    x + y + z
}

pub fn vector_dot(x1: f32, y1: f32, z1: f32, x2: f32, y2: f32, z2: f32) -> f32 {
    vector_sum_(vector_mul(x1, y1, z1, x2, y2, z2))
}

pub fn vector_dot_((x1, y1, z1): (f32, f32, f32), (x2, y2, z2): (f32, f32, f32)) -> f32 {
    vector_dot(x1, y1, z1, x2, y2, z2)
}

pub fn iter_integrate<'a>(
    time_step: f32,
    velocities: &'a [f32],
    positions: &'a [f32],
) -> impl Iterator<Item = f32> + 'a {
    positions
        .iter()
        .zip(velocities)
        .map(move |(position, velocity)| position + time_step * velocity)
}

#[cfg(test)]
mod tests {
    #[test]
    fn vector_mul() {
        assert_eq!(super::vector_mul(1., 2., 3., 0., 0., 0.), (0., 0., 0.));
        assert_eq!(super::vector_mul(1., 2., 3., 1., 1., 1.), (1., 2., 3.));
    }
    #[test]
    fn vector_dot() {
        assert_eq!(super::vector_dot(1., 2., 3., 1., 1., 1.), 6.);
        assert_eq!(super::vector_dot(1., 2., 3., 0., 0., 0.), 0.);
        assert_eq!(super::vector_dot(1., 2., 3., 1., 2., 3.), 14.);
        assert_eq!(super::vector_dot(1., 2., 3., 4., 5., 6.), 32.);
    }
    #[test]
    fn vector_dot_() {
        assert_eq!(super::vector_dot_((1., 2., 3.), (1., 1., 1.)), 6.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (0., 0., 0.)), 0.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (1., 2., 3.)), 14.);
        assert_eq!(super::vector_dot_((1., 2., 3.), (4., 5., 6.)), 32.);
    }
    #[test]
    fn vector_sub() {
        assert_eq!(super::vector_sub(1., 2., 3., 0., 0., 0.), (1., 2., 3.));
        assert_eq!(super::vector_sub(1., 2., 3., 1., 2., 3.), (0., 0., 0.));
    }
    #[test]
    fn vector_sub_() {
        assert_eq!(super::vector_sub_((1., 2., 3.), (0., 0., 0.)), (1., 2., 3.));
        assert_eq!(super::vector_sub_((1., 2., 3.), (1., 2., 3.)), (0., 0., 0.));
    }
    #[test]
    fn vector_add() {
        assert_eq!(super::vector_add(1., 2., 3., 0., 0., 0.), (1., 2., 3.));
        assert_eq!(super::vector_add(1., 2., 3., 4., 5., 6.), (5., 7., 9.));
    }
    #[test]
    fn vector_add_() {
        assert_eq!(super::vector_add_((1., 2., 3.), (0., 0., 0.)), (1., 2., 3.));
        assert_eq!(super::vector_add_((1., 2., 3.), (4., 5., 6.)), (5., 7., 9.));
    }
    #[test]
    fn vector_sum_() {
        assert_eq!(super::vector_sum_((1., 2., 3.)), 6.);
    }
    #[test]
    fn vector_scalar_mul() {
        assert_eq!(super::vector_scalar_mul(2., 1., 2., 3.), (2., 4., 6.));
    }
    #[test]
    fn iter_integrate() {
        let velocities = vec![10., 10.];
        let positions = vec![0., 15.];
        let new_positions: Vec<_> = super::iter_integrate(2., &velocities, &positions).collect();
        assert_eq!(new_positions, vec![20., 35.]);
    }
}
