//=============================================================================
//
//   Exercise code for the lecture
//   "Advanced Computer Graphics"
//
//   Adapted from Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2013 LGG, epfl
//
//   DO NOT REDISTRIBUTE
//=============================================================================


//== INCLUDES =================================================================

#include "Mass_spring_viewer.h"
#include "utils/gl_text.h"
#include <sstream>


#include <fstream>


//== IMPLEMENTATION ==========================================================


Mass_spring_viewer::
Mass_spring_viewer(const char* _title, int _width, int _height)
: Viewer_2D(_title, _width, _height)
{
    integration_         = Euler;
    collisions_          = Force_based;
    external_force_      = None;
    animate_             = false;
    area_forces_         = false;
    equilibrium_forces_  = false;
    log_kinetic_energy   = false;
    show_forces_         = false;

    time_step_           = 0.001;
    particle_radius_     = 0.03;

    particle_mass_       = 0.1;
    damping_             = 0.1;
    collision_stiffness_ = 1000.0;
    spring_stiffness_    = 1000.0;
    spring_damping_      = 1.0;
    area_stiffness_      = 100000.0;

    mouse_spring_.active = false;
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::keyboard(int key, int x, int y)
{
    switch (key)
    {
        // setup problem 1
        case '1':
        {
            body_.clear();
            body_.add_particle( vec2(-0.5, -0.5), vec2(14.0, -2.0), particle_mass_, false );
            glutPostRedisplay();
            break;
        }

        // setup problem 2
        case '2':
        {
            body_.clear();
            for (int i=0; i<100; ++i)
            {
                body_.add_particle( vec2(0.9* cos(i/50.0*M_PI), 0.9*sin(i/50.0*M_PI)), vec2(-sin(i/50.0*M_PI), cos(i/50.0*M_PI)), particle_mass_, false );
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 3
        case '3':
        {
            body_.clear();

            for (int i=0; i<10; ++i)
            {
                body_.add_particle( vec2(i*0.1, 0.8), vec2(0.0, 0.0), particle_mass_, i==0 );
            }

            for (unsigned int i=0; i<9; ++i)
            {
                body_.add_spring(i, i+1);
            }

            glutPostRedisplay();
            break;
        }

        // setup problem 4
        case '4':
        {
            body_.clear();

            body_.add_particle( vec2(-0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.0, 0.6), vec2(0.0, 0.0), particle_mass_, false );
            body_.add_particle( vec2( 0.1, 0.7), vec2(0.0, 0.0), particle_mass_, false );

            body_.add_spring(0, 1);
            body_.add_spring(0, 2);
            body_.add_spring(1, 2);

            body_.add_triangle(0, 1, 2);

            glutPostRedisplay();
            break;
        }

        // setup problem 5
        case '5':
        {
            body_.clear();

            for (int i=0; i<8; ++i)
            {
                body_.add_particle( vec2(-0.5+0.2*cos(0.25*i*M_PI), -0.5+0.2*sin(0.25*i*M_PI)), vec2(5.0, 5.0), particle_mass_, false );
            }

            body_.add_particle( vec2(-0.5, -0.5), vec2(5.0, 5.0), particle_mass_, false );

            for (unsigned int i=0; i<8; ++i)
            {
                body_.add_spring( i, (i+1) % 8 );
                body_.add_spring(i, 8);
                body_.add_triangle(i, (i+1)%8, 8);
            }

            glutPostRedisplay();
            break;
        }

        // switch between time integration methods
        case 'i':
        {
            switch (integration_)
            {
                case Euler:    integration_ = Midpoint; break;
                case Midpoint: integration_ = Verlet;   break;
                case Verlet:   integration_ = Implicit;    break;
                case Implicit: integration_ = Euler; break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between center and gravitation force
        case 'f':
        {
            switch (external_force_)
            {
                case None:        external_force_ = Center; break;
                case Center:      external_force_ = Gravitation; break;
                case Gravitation: external_force_ = None;      break;
            }
            glutPostRedisplay();
            break;
        }

        // switch between force-based and impulse-based collisions
        case 'c':
        {
            switch (collisions_)
            {
                case Force_based:   collisions_ = Impulse_based; break;
                case Impulse_based: collisions_ = Force_based;   break;
            }
            glutPostRedisplay();
            break;
        }


        // toggle area forces on/off
        case 'a':
        {
            area_forces_ = !area_forces_;
            glutPostRedisplay();
            break;
        }


        // visualization of particle forces on/off
        case 'v':
        {
            show_forces_ = !show_forces_;
            glutPostRedisplay();
            break;
        }

        case 'e':
        {
            equilibrium_forces_ = !equilibrium_forces_;
            glutPostRedisplay();
            break;
        }


        case 'l':
        {
            log_kinetic_energy = !log_kinetic_energy;

            if (log_kinetic_energy) {
                const char * name;
                switch (integration_)
                {
                    case Euler:    name = "Euler.txt"; break;
                    case Midpoint: name = "Midpoint.txt";   break;
                    case Verlet:   name = "Verlet.txt";    break;
                    case Implicit: name = "Implicit.txt"; break;
                }
                log_file.open(name, std::fstream::out);
            } else {
                log_file.close();
            }

            glutPostRedisplay();
            break;
        }


        // let parent class do the work
        default:
        {
            Viewer_2D::keyboard(key, x, y);
            break;
        }
    }
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::draw()
{
    // parent's status text
    Viewer_2D::draw();


    // draw some status text
    glDisable(GL_LIGHTING);
    glColor3f(1,0,0);
    std::ostringstream oss;

    oss.str("");
    oss << "Integration: ";
    switch (integration_)
    {
        case Euler:    oss << "Euler";    break;
        case Midpoint: oss << "Midpoint"; break;
        case Verlet:   oss << "Verlet";   break;
        case Implicit: oss << "Implicit"; break;
    }
    glText(20, height_-40, oss.str());

    oss.str("");
    oss << "#Particles: " << body_.particles.size();
    glText(20, height_-60, oss.str());

    oss.str("");
    oss << "#Springs: " << body_.springs.size();
    glText(20, height_-80, oss.str());

    oss.str("");
    oss << "#Triangles: " << body_.triangles.size();
    glText(20, height_-100, oss.str());

    oss.str("");
    oss << "Area Forces: " << (area_forces_ ? "on" : "off");
    glText(20, height_-120, oss.str());

    oss.str("");
    oss << "Collisions: " << (collisions_ == Force_based ? "force" : "impulse");
    glText(20, height_-140, oss.str());

    oss.str("");
    oss << "External force: ";
    switch (external_force_)
    {
        case None:        oss << "None";        break;
        case Center:      oss << "Center";      break;
        case Gravitation: oss << "Gravitation"; break;
    }
    glText(20, height_-160, oss.str());

    oss.str("");
    oss << "Visualize forces: " << (show_forces_ ? "on" : "off");
    glText(20, height_-180, oss.str());

    oss.str("");
    oss << "Equilibrium Forces: " << (equilibrium_forces_ ? "on" : "off");
    glText(20, height_-200, oss.str());

    oss.str("");
    oss << "Log kinetic energy: " << (log_kinetic_energy ? "on" : "off");
    glText(20, height_-220, oss.str());

    // draw walls
    glDisable(GL_LIGHTING);
    glLineWidth(1.0);
    glColor3f(0.5,0.5,0.5);
    glBegin(GL_LINE_STRIP);
    glVertex2f( -1.0,  1.0 );
    glVertex2f( -1.0, -1.0 );
    glVertex2f(  1.0, -1.0 );
    glVertex2f(  1.0,  1.0 );
    glVertex2f( -1.0,  1.0 );
    glEnd();


    // draw mouse spring
    if (mouse_spring_.active)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(5.0);
        glColor3f(1,0,0);
        glBegin(GL_LINES);
        glVertex2fv( body_.particles[mouse_spring_.particle_index].position.data() );
        glVertex2fv( mouse_spring_.mouse_position.data() );
        glEnd();
    }


    // draw particles, springs, triangles
    body_.draw(particle_radius_, show_forces_);
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::mouse(int _button, int _state, int _x, int _y)
{
    // need particles to do interaction
    if (!body_.particles.empty())
    {
        // mouse button release destroys current mouse spring
        if (_state == GLUT_UP)
        {
            mouse_spring_.active = false;
        }

        // mouse button press generates new mouse spring
        else if (_state == GLUT_DOWN)
        {
            // get point under mouse cursor
            vec2 p = pick(_x, _y);

            // find closest particle
            int   pidx = -1;
            float dmin = FLT_MAX;
            for (unsigned int i=0; i<body_.particles.size(); ++i)
            {
                float d = norm(p - body_.particles[i].position);
                if (d < dmin)
                {
                    dmin = d;
                    pidx = i;
                }
            }

            // construct mouse spring
            mouse_spring_.mouse_position = p;
            mouse_spring_.particle_index = pidx;
            mouse_spring_.active = true;
        }
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::motion(int _x, int _y)
{
    if (mouse_spring_.active)
    {
        // update mouse positions
        mouse_spring_.mouse_position = pick(_x, _y);
    }

    glutPostRedisplay();
}


//-----------------------------------------------------------------------------


void Mass_spring_viewer::time_integration(float dt)
{
    switch (integration_)
    {
        case Euler:
        {
            /** \todo (Part 1) Implement Euler integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */
            compute_forces ();
            for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
              if (!p_it->locked)
              {
                /// 1)
                p_it->acceleration = p_it->force / p_it->mass;

                /// 2)
                p_it->position = p_it->position + dt * p_it->velocity;

                /// 3)
                p_it->velocity = p_it->velocity + dt * p_it->acceleration;
              }

            break;
        }

        case Midpoint:
        {
            /** \todo (Part 2) Implement the Midpoint time integration scheme
             \li The Particle class has variables position_t and velocity_t to store current values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */
             compute_forces();
            for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
              if (!p_it->locked)
              {
                p_it->position_t = p_it->position;
                p_it->velocity_t = p_it->velocity;
                p_it->acceleration = p_it->force/p_it->mass;

                p_it->position += 0.5*dt * p_it->velocity;
                p_it->velocity += 0.5*dt * p_it->acceleration;
              }
            compute_forces();
            for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
              if (!p_it->locked)
              {
                p_it->acceleration = p_it->force / p_it->mass;


                p_it->position = p_it->position_t + dt * p_it->velocity;
                p_it->velocity = p_it->velocity_t + dt * p_it->acceleration;
              }

            break;
        }


        case Verlet:
        {
            /** \todo (Part 2) Implement the Verlet time integration scheme
             \li The Particle class has a variable acceleration to remember the previous values
             \li Hint: compute_forces() computes all forces for the current positions and velocities.
             */            
            compute_forces();
            // x(t)  v(t-h)
            for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
              if (!p_it->locked)
              {
                vec2 a = p_it->force / p_it->mass;
                p_it->velocity += 0.5*dt*(p_it->acceleration + a);

                // x(t) v(t)
                p_it->position += dt*p_it->velocity + 0.5*dt*dt*a;
                p_it->acceleration = a;
              }

            break;
        }

        case Implicit:
        {
            /// The usual force computation method is called, and then the jacobian matrix dF/dx is calculated
            compute_forces ();
            compute_jacobians ();

            /// Finally the linear system is composed and solved
            solver_.solve (dt, particle_mass_, body_.particles);

            break;
        }
    }


    // impulse-based collision handling
    if (collisions_ == Impulse_based)
    {
        impulse_based_collisions();
    }


    // E_c = 1/2 m*v^2
    float total_E = 0.0f;
    for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        total_E += 0.5 * p_it->mass * dot(p_it->velocity, p_it->velocity);


    log_file << total_E << "\n";
    glutPostRedisplay();
}


//-----------------------------------------------------------------------------

void
Mass_spring_viewer::compute_forces()
{
    // clear forces
    for (unsigned int i=0; i<body_.particles.size(); ++i)
        body_.particles[i].force = vec2(0,0);


    /** \todo (Part 1) Implement center force
     */
    if (external_force_ == Center)
    {
      float c = 10.f;
      for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        p_it->force += c * (vec2 (0.f, 0.f) - p_it->position);
    }


    /** \todo (Part 1) Implement damping force
     \li The damping coefficient is given as member damping_
     */
    /// Do damping only for explicit methods
    if (integration_ != Implicit)
      for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        p_it->force -= damping_ * p_it->velocity;


    /** \todo (Part 1) Implement gravitation force
     \li Particle mass available as particle_mass_
     */
    if (external_force_ == Gravitation)
    {
      for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        p_it->force += 9.81 * p_it->mass * vec2 (0.f, -1.f);
    }


    /** \todo (Part 1) Implement force based boundary collisions
     \li Collision coefficient given as collision_stiffness_
     */
    // collision forces
    if (collisions_ == Force_based)
    {
        float planes[4][3] = {
            {  0.0,  1.0, 1.0 },
            {  0.0, -1.0, 1.0 },
            {  1.0,  0.0, 1.0 },
            { -1.0,  0.0, 1.0 }
        };

        for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        {
          for (size_t plane_i = 0; plane_i < 4; ++plane_i)
          {
            float d = planes[plane_i][0] * p_it->position[0] +
                      planes[plane_i][1] * p_it->position[1] +
                      planes[plane_i][2] - particle_radius_;

            if (d < 0.f)
              p_it->force += collision_stiffness_ * (-d) * vec2 (planes[plane_i][0],
                                                                 planes[plane_i][1]);
          }
        }
    }


    /** \todo (Part 1) Compute force of the interactive mass spring
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */
    if (mouse_spring_.active)
    {
        Particle& particle = body_.particles[ mouse_spring_.particle_index ];

        vec2 pos0 = particle.position;
        vec2 pos1 = mouse_spring_.mouse_position;

        particle.force += -(spring_stiffness_ * norm (pos0 - pos1) +
                            spring_damping_ * dot (particle.velocity, pos0 - pos1) / norm (pos0 - pos1)) * (pos0 - pos1) / norm (pos0 - pos1);
    }


    /** \todo (Part 1) Compute spring forces
     \li Required information about springs in the scene are found in body_.springs
     \li Required coefficients are given as spring_stiffness_ and spring_damping_
     */

    for (std::vector<Spring>::iterator s_it = body_.springs.begin (); s_it != body_.springs.end (); ++s_it)
    {
      const vec2 &pos0 = s_it->particle0->position;
      const vec2 &pos1 = s_it->particle1->position;
      const vec2 &vel0 = s_it->particle0->velocity;
      const vec2 &vel1 = s_it->particle1->velocity;

      vec2 force = -(spring_stiffness_ * (norm (pos0 - pos1) - s_it->rest_length) +
                     spring_damping_ * dot (vel0 - vel1, pos0 - pos1) / norm (pos0 - pos1)) * (pos0 - pos1) / norm (pos0 - pos1);
      s_it->particle0->force += force;
      s_it->particle1->force += -force;
    }



    /** \todo (Part 2) Compute more forces in part 2 of the exercise: triangle-area forces, binding forces, etc.
     */
    if (area_forces_) {
        for (std::vector<Triangle>::iterator t_it = body_.triangles.begin (); t_it != body_.triangles.end (); ++t_it)
        {
            /*
            F_i = - k_A * (area - A) * d_area/d_i


            */

            const vec2& p0 = t_it->particle0->position;
            const vec2& p1 = t_it->particle1->position;
            const vec2& p2 = t_it->particle2->position;
            vec2 d_area_0 = perp(0.5 * (p1 - p2));
            vec2 d_area_1 = perp(0.5 * (p2 - p0));
            vec2 d_area_2 = perp(0.5 * (p0 - p1));

            float tmp = - area_stiffness_ * (t_it->area() - t_it->rest_area);

            t_it->particle0->force += tmp * d_area_0;
            t_it->particle1->force += tmp * d_area_1;
            t_it->particle2->force += tmp * d_area_2;
        }
    }


    if (equilibrium_forces_) {
        for (std::vector<Particle>::iterator p_it1 = body_.particles.begin (); p_it1 != body_.particles.end (); ++p_it1) {
            for (std::vector<Particle>::iterator p_it2 = p_it1+1; p_it2 != body_.particles.end (); ++p_it2) {
                const vec2 &p1 = p_it1->position;
                const vec2 &p2 = p_it2->position;
                float dist = distance(p1, p2);
                vec2 d = p1 - p2;
                float f = 0.1f / (dist*dist);

                p_it1->force += f*d;
                p_it2->force += -f*d;
            }
        }

    }
}




//-----------------------------------------------------------------------------


void Mass_spring_viewer::impulse_based_collisions()
{
    /** \todo (Part 2) Handle collisions based on impulses
     */
     float planes[4][3] = {
            {  0.0,  1.0, 1.0 },
            {  0.0, -1.0, 1.0 },
            {  1.0,  0.0, 1.0 },
            { -1.0,  0.0, 1.0 }
        };

        for (std::vector<Particle>::iterator p_it = body_.particles.begin (); p_it != body_.particles.end (); ++p_it)
        {
          for (size_t plane_i = 0; plane_i < 4; ++plane_i)
          {
            float d = planes[plane_i][0] * p_it->position[0] +
                      planes[plane_i][1] * p_it->position[1] +
                      planes[plane_i][2] - particle_radius_;

            if (d < 0.f)
              p_it->velocity += collision_stiffness_ * (-d) * vec2 (planes[plane_i][0],
                                                                 planes[plane_i][1]);
          }
        }
}
//=============================================================================

void Mass_spring_viewer::compute_jacobians ()
{
  /// Clear the solver matrices
  solver_.clear ();

  /** \todo (Part 2) Implement the corresponding jacobians for each of the force types.
   * Use the code from compute_forces() as the starting ground.
   */
}

//=============================================================================
void ImplicitSolver::solve (float dt, float mass,
                            std::vector<Particle> &particles)
{
  int num_particles = particles.size ();

  /// Build the Jacobian matrix from the sparse set of elements
  Eigen::SparseMatrix<float> J (2 * num_particles, 2 * num_particles);
  J.setFromTriplets (triplets_.begin (), triplets_.end ());

  /// Build up the position, velocity and force vectors
  Eigen::VectorXf pos_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  velocity_vec (Eigen::VectorXf::Zero (2 * num_particles)),
                  force_vec (Eigen::VectorXf::Zero (2 * num_particles));
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    pos_vec (2 * p_i + 0) = particles[p_i].position[0];
    pos_vec (2 * p_i + 1) = particles[p_i].position[1];
    velocity_vec (2 * p_i + 0) = particles[p_i].velocity[0];
    velocity_vec (2 * p_i + 1) = particles[p_i].velocity[1];
    force_vec (2 * p_i + 0) = particles[p_i].force[0];
    force_vec (2 * p_i + 1) = particles[p_i].force[1];
  }

  /// Kick out the fixed particles by creating a sparse selection matrix
  std::vector<Eigen::Triplet<float> > triplets_selection;
  int valid_particle_index = 0;
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
    if (!particles[p_i].locked)
    {
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 0, 2 * valid_particle_index + 0, 1.f));
      triplets_selection.push_back (Eigen::Triplet<float> (2 * p_i + 1, 2 * valid_particle_index + 1, 1.f));
      valid_particle_index ++;
    }
  Eigen::SparseMatrix<float> mat_selection (2 * num_particles, 2 * valid_particle_index);
  mat_selection.setFromTriplets (triplets_selection.begin (), triplets_selection.end ());

  /// Sparse identity matrix
  Eigen::SparseMatrix<float> Id (2 * valid_particle_index, 2 * valid_particle_index);
  Id.setIdentity ();

  /// Apply the selection matrix on each vector and the Jacobian matrix
  pos_vec = mat_selection.transpose () * pos_vec;
  velocity_vec = mat_selection.transpose () * velocity_vec;
  force_vec = mat_selection.transpose () * force_vec;
  J = mat_selection.transpose () * J * mat_selection;

  /// Build the right and left hand sides of the linear system
  Eigen::SparseMatrix<float> A = Id - dt * dt / mass * J;
  Eigen::VectorXf b;
  b = dt * velocity_vec + dt * dt / mass * force_vec + (Id - dt * dt / mass * J) * pos_vec;

  /// Solve the system and use the selection matrix again to arrange the new positions in a vector
  linear_solver_.analyzePattern (A);
  linear_solver_.compute (A);
  Eigen::VectorXf new_pos = mat_selection * linear_solver_.solve (b);

  /// Extract the positions from the solution vector and set the new positions and velocities inside the particle structures
  for (size_t p_i = 0; p_i < num_particles; ++p_i)
  {
    if (!particles[p_i].locked)
    {
      vec2 pos_update (new_pos (2 * p_i + 0), new_pos (2 * p_i + 1));
      particles[p_i].velocity = (pos_update - particles[p_i].position) / dt;
      particles[p_i].position = pos_update;
    }
  }
}
