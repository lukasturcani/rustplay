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

pub struct PhysicsConfig {
    pub time_step: f64,
    pub max_x: f64,
    pub max_y: f64,
    pub max_z: f64,
}
