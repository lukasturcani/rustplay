use kiss3d::nalgebra::Point3;
use kiss3d::window::Window;

pub fn draw_box(color: &Point3<f32>, x: f32, y: f32, z: f32, window: &mut Window) {
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(x, 0., 0.), color);
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(0., y, 0.), color);
    window.draw_line(&Point3::new(0., 0., 0.), &Point3::new(0., 0., z), color);
    window.draw_line(&Point3::new(x, 0., 0.), &Point3::new(x, y, 0.), color);
    window.draw_line(&Point3::new(x, 0., 0.), &Point3::new(x, 0., z), color);
    window.draw_line(&Point3::new(x, y, 0.), &Point3::new(0., y, 0.), color);
    window.draw_line(&Point3::new(x, y, 0.), &Point3::new(x, y, z), color);
    window.draw_line(&Point3::new(0., y, 0.), &Point3::new(0., y, z), color);
    window.draw_line(&Point3::new(0., 0., z), &Point3::new(x, 0., z), color);
    window.draw_line(&Point3::new(0., 0., z), &Point3::new(0., y, z), color);
    window.draw_line(&Point3::new(x, 0., z), &Point3::new(x, y, z), color);
    window.draw_line(&Point3::new(x, y, z), &Point3::new(0., y, z), color);
}
