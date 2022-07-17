pub struct PhysicsConfig {}

#[derive(Debug)]
pub struct Quaternion {}

#[derive(Debug)]
pub struct PhysicsData {
    pub positions_x: Vec<f32>,
    pub positions_y: Vec<f32>,
    pub positions_z: Vec<f32>,
    pub velocities_x: Vec<f32>,
    pub velocities_y: Vec<f32>,
    pub velocities_z: Vec<f32>,
    pub orientations: Vec<Quaternion>,
}
