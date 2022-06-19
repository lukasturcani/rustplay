use std::vec::Vec;

pub mod vec3;

pub fn loop_integrate(time_step: f64, positions: &[f64], velocities: &[f64]) -> Vec<f64> {
    let mut output = Vec::with_capacity(positions.len());
    for i in 0..positions.len() {
        output.push(positions[i] + time_step * velocities[i]);
    }
    output
}

pub fn iter_integrate<'a>(
    time_step: f64,
    velocities: &'a [f64],
    positions: &'a [f64],
) -> impl Iterator<Item = f64> + 'a {
    let time_step_ = time_step.clone();
    positions
        .iter()
        .zip(velocities)
        .map(move |(position, velocity)| position + time_step_ * velocity)
}

pub fn loop_apply_bounds<'a>(
    min: f64,
    max: f64,
    velocities: &mut [f64],
    positions: &'a mut [f64],
) -> &'a mut [f64] {
    for i in 0..positions.len() {
        if positions[i] < min {
            positions[i] = min + min - positions[i];
            velocities[i] *= -1.0;
        } else if positions[i] > max {
            positions[i] = max + max - positions[i];
            velocities[i] *= -1.0;
        }
    }
    positions
}


pub fn iter_apply_bounds<'a>(
    min: f64,
    max: f64,
    velocities: &mut [f64],
    positions: &'a mut [f64],
) -> &'a mut [f64] {

    positions.iter_mut().zip(velocities).for_each(
        |(position, velocity)| {
            if *position < min {
                *position = &min + &min - *position;
                *velocity *= -1.0;
            } else if *position > max {
                *position = max + max - *position;
                *velocity *= -1.0;
            }
        }
    );
    positions
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
