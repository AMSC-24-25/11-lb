#include "Node.hpp"
#include "Lattice.hpp"
#include "Auxiliary.cpp"

Node::Node(const std::vector<int>& boundary_, const bool obs_, const std::vector<int>& position_)
            : boundary(boundary_), obstacle(obs_), position(position_)
{
    velocity.resize(2, 0.0);
    f.resize(9,0.0);
    newF = f;  
}

void Node::updateMacroscopic()
{
    //if the node is an obstacle nothing to do 
    if (obstacle == true) {return;}
    else
    {
        density = 0.0;
        for (int i = 0; i < f.size(); ++i)
        {
            density += f[i];
        }
        for (int j=0; j<2; j++)
        {
            velocity.at(j) = 0.0;
            for (int k=0; k<9; k++)
            {
                velocity.at(j) += f.at(k)*dir_by_dim.at(j).at(k);
            }
            velocity.at(j)/=density;
        }
    }
}

void Node::equilibriumCollision(double omegaP, double halfOmegaSum, double halfOmegaSub)
{
    if (obstacle == true) {return;}
    else 
    {
        // equilibrium
        const double vel_sq = dotProduct(velocity, velocity);
        const double k1 = 3/2*vel_sq;
        std::vector<double> f_eq(9,0.0);
        for (int i=0 ; i<2 ; i++)
        {
            const double k2 = 3*dotProduct(dir.at(i), velocity);
            f_eq.at(i) = density*w.at(i)*(1+k2+1/2*k2*k2-k1);
        }
        // Collision, index 0 calculated outside the cicle
        newF.at(0) = omegaP*f_eq.at(0)+f.at(0)*(1-omegaP);
        for (int i=1; i<9; i++)
        {
            newF.at(i) = halfOmegaSum*f_eq.at(i) + halfOmegaSub*f_eq.at(opp.at(i)) 
                            + (1.0-halfOmegaSum)*f.at(i) - halfOmegaSub*f.at(opp.at(i));
        }
    }
}

void Node::initializeEquilibrium()
{
    if (obstacle==true)
    {
        velocity.at(0) = 0.0;
        velocity.at(1) = 0.0;
        return;
    }
    else 
    {
        const double vel_sq = dotProduct(velocity, velocity);
        const double k1 = 3/2*vel_sq;
        std::vector<double> f_eq(9,0.0);
        for (int i=0 ; i<2 ; i++)
        {
            const double k2 = 3*dotProduct(dir.at(i), velocity);
            f.at(i) = density*w.at(i)*(1+k2+1/2*k2*k2-k1);
        }
    }
}

void Node::streaming( Lattice& lattice )
{
    if (obstacle == true) {return;}
    else 
    {
        // Stream for index 0 and then for other indices in the for loop 
        f.at(0) = newF.at(0);
        std::vector<int> newPos(2);
        const unsigned int NX = lattice.node_matrix.shape().at(0);
        const unsigned int NY = lattice.node_matrix.shape().at(1);
        for (int i=1; i<9; i++)
        {
            // Evaluate new position
            for(int j=0; j<2; j++)
            {
                newPos.at(j) = position.at(j) + (int)dir.at(i).at(j);
            }
            // If not an obstacle and is inside lattice compute the streaming
            if ( newPos.at(0)>=0 && newPos.at(0)<NX      // in rows
                && newPos.at(1)>=0 && newPos.at(1)<NY    // in columns
                && lattice.node_matrix(newPos.at(0),newPos.at(1)).isObstacle()==false ) // is fluid  
            {
                lattice.node_matrix(newPos.at(0),newPos.at(1)).setFAtIndex(i,newF.at(i));
            }
        }
    }
}

void Node::applyInletBoundary(Lattice& lattice, std::vector<double> boundary_velocity)
{
    if (obstacle == true) {return;}
    else if (boundary_velocity.size()!=2) {return;}
    else 
    {
        const unsigned int NX = lattice.node_matrix.shape().at(0);
        const unsigned int NY = lattice.node_matrix.shape().at(1);
        // Default zero velocity on walls
        if ( position.at(0)==0 || position.at(0)==(NX-1) || position.at(1)==0 || position.at(1)==(NY-1))
        {
            velocity.at(0) = 0.0;
            velocity.at(1) = 0.0;
        } 
        // WindTunnel-Like problem left
        if (position.at(0)==0)
        {
            velocity.at(0) = boundary_velocity.at(0);
            velocity.at(1) = boundary_velocity.at(1);
        }
    }
}

void Node::applyZouHeBoundary(Lattice &lattice)
{
    if (obstacle == true) {return;}
    else 
    {
        const unsigned int NX = lattice.node_matrix.shape().at(0);
        const unsigned int NY = lattice.node_matrix.shape().at(1);

        // BC on top wall
        if ( position.at(0)!=0 && position.at(0)!=(NX-1) && position.at(1)==0 )
        {
            density = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(2) + f.at(5) + f.at(6))) / (1.0 + velocity.at(1));
            f.at(4) = f.at(2) - 2.0 / 3.0 * density * velocity.at(1);
            f.at(7) = f.at(5) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * density * velocity.at(0) - 1.0 / 6.0 * density * velocity.at(1);
            f.at(8) = f.at(6) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * density * velocity.at(0) - 1.0 / 6.0 * density * velocity.at(1);
        }
        // BC on right wall
        else if (position.at(0)==(NX-1) && position.at(1)!=0 && position.at(1)!=(NY-1))
        {
            // I know rho so I evaluate velocity 
            density = 1;
            velocity.at(0) = f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(1) + f.at(5) + f.at(8)) - 1.0;
            f.at(3) = f.at(1) - 2.0 / 3.0 * velocity.at(0);
            f.at(6) = f.at(8) - 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * velocity.at(0);
            f.at(7) = f.at(5) + 0.5 * (f.at(2) - f.at(4)) - 1.0 / 6.0 * velocity.at(0);
        }
        // BC on bottom wall
        else if (position.at(0)!=0 && position.at(0)!=(NX-1) && position.at(1)==(NY-1))
        {
            density = (f.at(0) + f.at(1) + f.at(3) + 2.0 * (f.at(4) + f.at(7) + f.at(8))) / (1.0 - velocity.at(1));
            f.at(2) = f.at(4) + 2.0 / 3.0 * density * velocity.at(1);
            f.at(5) = f.at(7) - 0.5 * (f.at(1) - f.at(3)) + 0.5 * density * velocity.at(0) + 1.0 / 6.0 * density * velocity.at(1);
            f.at(6) = f.at(8) + 0.5 * (f.at(1) - f.at(3)) - 0.5 * density * velocity.at(0) + 1.0 / 6.0 * density * velocity.at(1);
        }
        // BC on left wall
        else if (position.at(0)==0 && position.at(1)!=0 && position.at(1)!=(NY-1))
        {
            density = (f.at(0) + f.at(2) + f.at(4) + 2.0 * (f.at(3) + f.at(6) + f.at(7))) / (1.0 - velocity.at(0));
            f.at(1) = f.at(3) + 2.0 / 3.0 *density* velocity.at(0);
            f.at(5) = f.at(7) - 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * density * velocity.at(0) + 0.5 * density * velocity.at(1);
            f.at(8) = f.at(6) + 0.5 * (f.at(2) - f.at(4)) + 1.0 / 6.0 * density * velocity.at(0) - 0.5 * density * velocity.at(1);
        }
        // BC on top right corner
        else if (position.at(0)==(NX-1) && position.at(1)==0)
        {
            density = lattice.node_matrix(position.at(0)-1,position.at(1)).getDensity();

            f.at(3) = f.at(1) - 2.0 / 3.0 * density * velocity.at(0);
            f.at(4) = f.at(2) - 2.0 / 3.0 * density * velocity.at(1);
            f.at(7) = f.at(5) - 1.0 / 6.0 * density * velocity.at(0) - 1.0 / 6.0 * density * velocity.at(1);
            f.at(8) = 0;
            f.at(6) = 0;
            f.at(0) = density - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
        }
        // BC on bottom right corner
        else if (position.at(0)==(NX-1) && position.at(1)==(NY-1))
        {
            density = lattice.node_matrix(position.at(0)-1, position.at(1)).getDensity();
    
            f.at(3) = f.at(1) - 2.0 / 3.0 * density * velocity.at(0);
            f.at(2) = f.at(4) + 2.0 / 3.0 * density * velocity.at(1);
            f.at(6) = f.at(8) + 1.0 / 6.0 * density * velocity.at(1) - 1.0 / 6.0 * density * velocity.at(0);
            f.at(7) = 0;
            f.at(5) = 0;
            f.at(0) = density - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
        }
        // BC on bottom left corner
        else if (position.at(0)==0 && position.at(1)==(NY-1))
        {
            density = lattice.node_matrix(position.at(0)+1, position.at(1)).getDensity();
    
            f.at(1) = f.at(3) + 2.0 / 3.0 * density * velocity.at(0);
            f.at(2) = f.at(4) + 2.0 / 3.0 * density * velocity.at(1);
            f.at(5) = f.at(7) + 1.0 / 6.0 * density * velocity.at(0) + 1.0 / 6.0 * density * velocity.at(1);
            f.at(6) = 0;
            f.at(8) = 0;
            f.at(0) = density - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(5) - f.at(7);
        }
        // BC on top left corner
        else if (position.at(0)==0 && position.at(1)==0)
        {
            density = lattice.node_matrix(position.at(0)+1, position.at(1)).getDensity();

            f.at(1) = f.at(3) + 2.0 / 3.0 * density * velocity.at(0);
            f.at(4) = f.at(2) - 2.0 / 3.0 * density * velocity.at(1);
            f.at(8) = f.at(6) + 1.0 / 6.0 * density * velocity.at(0) - 1.0 / 6.0 * density * velocity.at(1);
            f.at(7) = 0;
            f.at(5) = 0;
            f.at(0) = density - f.at(1) - f.at(2) - f.at(3) - f.at(4) - f.at(6) - f.at(8);
        }
    }
}

void Node::bounceBack()
{
    if (obstacle == true) {return;}
    else 
    {
        // Check of the boundary position wrt the actual Node
        if (boundary.at(0)==1) f.at(3) = newF.at(1);
        if (boundary.at(0)==-1) f.at(1) = newF.at(3);
        if (boundary.at(1)==1) f.at(2) = newF.at(4);
        if (boundary.at(1)==-1) f.at(4) = newF.at(2);
        if (boundary.at(2)==1) f.at(6) = newF.at(8);
        if (boundary.at(2)==-1) f.at(8) = newF.at(6);
        if (boundary.at(3)==1) f.at(5) = newF.at(7);
        if (boundary.at(3)==-1) f.at(7) = newF.at(5);
    }
}

void Node::computeDragAndLift(double &Cx, double &Cy, double rhoRef, double LRef, double URef)
{
    if (obstacle || (boundary.at(0) == 0 && boundary.at(1) == 0 && boundary.at(2) == 0 && boundary.at(3) == 0))
        return;

    double Fx = 0.0;
    double Fy = 0.0;

    if (boundary.at(0) == 1)
    {
        Fx += newF.at(1) + f.at(3);
    }
    else if (boundary.at(0) == -1)
    {
        Fx -= newF.at(3) + f.at(1);
    }
    if (boundary.at(1) == 1)
    {
        Fy += newF.at(4) + f.at(2);
    }
    else if (boundary.at(1) == -1)
    {
        Fy -= newF.at(2) + f.at(4);
    }
    if (boundary.at(2) == 1)
    {
        double t = newF.at(8) + f.at(6);
        Fx += t;
        Fy += t;
    }
    else if (boundary.at(2) == -1)
    {
        double t = newF.at(6) + f.at(8);
        Fx -= t;
        Fy -= t;
    }
    if (boundary.at(3) == 1)
    {
        double t = newF.at(7) + f.at(5);
        Fx -= t;
        Fy += t;
    }
    else if (boundary.at(3) == -1)
    {
        double t = newF.at(5) + f.at(7);
        Fx += t;
        Fy -= t;
    }

    // Normalize to get dimensionless drag/lift coefficients
    double normFactor = rhoRef * LRef * URef * URef;   // Sarebbe rho*v^2*S solo che in 2D quindi la superficie è una lunghezza 
    double localCx = -2.0 * Fx / normFactor;
    double localCy = -2.0
     * Fy / normFactor;

    Cx += localCx;
    Cy += localCy;
}

const double& Node::getDensity() const
{
    return density;
}

const std::vector<double>& Node::getVelocity() const
{
    return velocity;
}

const std::vector<double>& Node::getF( bool age) const
{
    // Age == true --> 'vecchia' f
    // Age == false  --> 'nuova' f
    if (age==true) return f;
    else return newF;
}

const std::vector<int>& Node::getBoundary() const
{
    return boundary;
}

bool Node::isObstacle() const
{
    return obstacle;
}

void Node::setFAtIndex( int index, const double& value )
{
    // Helper to set distribution function at a given value
    f.at(index) = value;
}
