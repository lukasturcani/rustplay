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
    fn iter_integrate() {
        let velocities = vec![10., 10.];
        let positions = vec![0., 15.];
        let new_positions: Vec<_> = super::iter_integrate(2., &velocities, &positions).collect();
        assert_eq!(new_positions, vec![20., 35.]);
    }
}
