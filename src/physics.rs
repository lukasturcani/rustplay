pub struct PhysicsConfig {
    pub time_step: f64,
    pub max_x: f64,
    pub max_y: f64,
    pub max_z: f64,
}

pub struct PhysicsData {
    pub positions_x: Vec<f64>,
    pub positions_y: Vec<f64>,
    pub positions_z: Vec<f64>,
    pub velocities_x: Vec<f64>,
    pub velocities_y: Vec<f64>,
    pub velocities_z: Vec<f64>,
}

pub fn get_random_physics_data(
    generator: &mut impl rand::Rng,
    num_objects: u64,
    max_position_x: f64,
    max_position_y: f64,
    max_position_z: f64,
    max_velocity_x: f64,
    max_velocity_y: f64,
    max_velocity_z: f64,
) -> PhysicsData {
    PhysicsData {
        positions_x: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_x)
            .collect(),
        positions_y: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_y)
            .collect(),
        positions_z: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_position_z)
            .collect(),
        velocities_x: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_x)
            .collect(),
        velocities_y: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_y)
            .collect(),
        velocities_z: (0..num_objects)
            .map(|_| generator.gen::<f64>() * max_velocity_z)
            .collect(),
    }
}

pub fn iter_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) {
    data.positions_x = iter_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    )
    .collect();
    iter_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );
    data.positions_y = iter_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    )
    .collect();
    iter_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    data.positions_z = iter_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    )
    .collect();
    iter_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
}

pub fn loop_take_time_step(config: &PhysicsConfig, data: &mut PhysicsData) {
    data.positions_x = loop_integrate(
        config.time_step,
        &data.velocities_x[..],
        &data.positions_x[..],
    );
    loop_apply_bounds(
        0.,
        config.max_x,
        &mut data.velocities_x[..],
        &mut data.positions_x[..],
    );
    data.positions_y = loop_integrate(
        config.time_step,
        &data.velocities_y[..],
        &data.positions_y[..],
    );
    loop_apply_bounds(
        0.,
        config.max_y,
        &mut data.velocities_y[..],
        &mut data.positions_y[..],
    );
    data.positions_z = loop_integrate(
        config.time_step,
        &data.velocities_z[..],
        &data.positions_z[..],
    );
    loop_apply_bounds(
        0.,
        config.max_z,
        &mut data.velocities_z[..],
        &mut data.positions_z[..],
    );
}

fn loop_integrate(time_step: f64, positions: &[f64], velocities: &[f64]) -> Vec<f64> {
    let mut output = Vec::with_capacity(positions.len());
    for i in 0..positions.len() {
        output.push(positions[i] + time_step * velocities[i]);
    }
    output
}

fn iter_integrate<'a>(
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

fn loop_apply_bounds<'a>(
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

fn iter_apply_bounds<'a>(
    min: f64,
    max: f64,
    velocities: &mut [f64],
    positions: &'a mut [f64],
) -> &'a mut [f64] {
    positions
        .iter_mut()
        .zip(velocities)
        .for_each(|(position, velocity)| {
            if *position < min {
                *position = &min + &min - *position;
                *velocity *= -1.0;
            } else if *position > max {
                *position = max + max - *position;
                *velocity *= -1.0;
            }
        });
    positions
}
