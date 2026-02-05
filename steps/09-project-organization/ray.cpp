#include "ray.h"
#include "constants.h"
#include <cmath>

// TODO: Implement Ray constructor
// - Initialize scenario and startFrame in member initializer list
// - Convert Cartesian (x, y) to polar (r, phi)
// - Store initialVelocityAngle from velocity direction
// - Convert Cartesian velocity to polar (v_r, v_phi)
// - Calculate conserved quantities (L, E)
// - Add initial position to trail

// TODO: Implement isCaptured() method
// - Return true if r <= BlackHole::rs * 1.01

// TODO: Implement hasEscaped() method
// - Return true if r > maxDistance

// TODO: Implement recordPosition() method
// - Convert polar (r, phi) to Cartesian (x, y)
// - Add position to trail

// TODO: Implement integrate() method
// - Check if ray should integrate (active, not captured, not escaped)
// - Call Physics::rk4Step (from physics module)
// - Record position
// - Update deflection

// TODO: Implement updateDeflection() method
// - Convert polar velocity to Cartesian velocity
// - Calculate current velocity direction angle
// - Compute deflection as change from initial velocity angle
// - Handle angle wrapping

// TODO: Implement isActive() method
// - Return true if currentFrame >= startFrame
