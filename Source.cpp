#include <iostream>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <random>

struct vec2
{
	float m_x, m_y;
};

template<typename T>
using Array = std::vector<T>;
using Vec2Array = Array<vec2>;

struct ParticleSystem
{
	Vec2Array m_aPositions;
	Vec2Array m_aVelocities;
	Array<float> m_aCharges;
	Array<float> m_aMasses;
	auto &addParticle(vec2 const &pos, float mass, float charge)
	{
		m_aPositions.push_back(pos);
		m_aVelocities.push_back({0.0f, 0.0f});
		m_aCharges.push_back(charge);
		m_aMasses.push_back(mass);
		return *this;
	}
};

struct Force
{
	float m_x = 0.0f, m_y = 0.0f;
};

Force getForce(vec2 const &pos1, float charge1, vec2 const &pos2, float charge2)
{
	float dx = pos1.m_x - pos2.m_x;
	float dy = pos1.m_y - pos2.m_y;
	float dist2 = dx * dx + dy * dy;
	static const float EPS = 1.0f;
	dist2 += EPS;
	float dist = sqrtf(dist2);
	float chargeMul = charge1 * charge2;
	float force = chargeMul / dist2;
	//here I make the force repulsive so that the particles wont collide
	float force_x = dx * (dist - 10.0f) / dist2 * force;
	float force_y = dy * (dist - 10.0f) / dist2 * force;
	return {force_x, force_y};
}

float getPotential(vec2 const &pos, vec2 const &parPos, float charge)
{
	float dx = parPos.m_x - pos.m_x;
	float dy = parPos.m_y - pos.m_y;
	float dist2 = dx * dx + dy * dy;
	static const float EPS = 1.0e-6f;
	float dist = sqrtf(dist2) + EPS;
	return charge / dist;
}

struct World
{
	int m_width, m_height;
	vec2 genRandomPoint() const
	{
		int x = rand() % m_width;
		int y = rand() % m_height;
		return {float(x), float(y)};
	}
	float genRandom() const
	{
		return float(rand()) / RAND_MAX;
	}
};

void iter( ParticleSystem &parSys, World const &world, float dt )
{
	struct Force
	{
		float m_x = 0.0f, m_y = 0.0f;
	};
	Array<Force> aForces;
	auto const &aPos = parSys.m_aPositions;
	aForces.resize(aPos.size());
	for (int i = 0; i < aPos.size(); i++)
	{
		auto par1 = aPos[i];
		for (int j = i + 1; j < aPos.size(); j++)
   		{
			auto par2 = aPos[j];
			auto force = getForce(par1, parSys.m_aCharges[i], par2, parSys.m_aCharges[j]);
			aForces[i].m_x += force.m_x;
			aForces[i].m_y += force.m_y;
			aForces[j].m_x += -force.m_x;
			aForces[j].m_y += -force.m_y;
		}
		// apply force from the imaginary particles. Just take all the particles and shift by the world size
		static int offsets[][2] = {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
		for (int offsetId = 0; offsetId < 8; offsetId++)
		{
			for (int j = 0; j < aPos.size(); j++)
	   		{
				auto par2 = aPos[j];
				par2.m_x += offsets[offsetId][0] * world.m_width;
				par2.m_y += offsets[offsetId][1] * world.m_height;
				auto force = getForce(par1, parSys.m_aCharges[i], par2, parSys.m_aCharges[j]);
				aForces[i].m_x += force.m_x;
				aForces[i].m_y += force.m_y;
			}
		}
	}
	for (int i = 0; i < aPos.size(); i++)
	{
		parSys.m_aVelocities[i].m_x += aForces[i].m_x / parSys.m_aMasses[i] * dt;
		parSys.m_aVelocities[i].m_y += aForces[i].m_y / parSys.m_aMasses[i] * dt;
		parSys.m_aPositions[i].m_x += parSys.m_aVelocities[i].m_x * dt;
		parSys.m_aPositions[i].m_y += parSys.m_aVelocities[i].m_y * dt;
	}
}

void regularize(ParticleSystem &parSys, World const &world)
{
	auto &aPos = parSys.m_aPositions;
	auto &aVel = parSys.m_aVelocities;
	for (int i = 0; i < aPos.size(); i++)
	{
		if (aPos[i].m_x < 0.0f)
		{
			aPos[i].m_x += world.m_width - 1;
		}
		if (aPos[i].m_x > world.m_width - 1)
		{
			aPos[i].m_x -= world.m_width - 1;
		}
		if (aPos[i].m_y < 0.0f)
	   	{
			aPos[i].m_y += world.m_height - 1;
		}
		if (aPos[i].m_y > world.m_height - 1)
   		{
			aPos[i].m_y -= world.m_height - 1;
		}
		aVel[i].m_x *= 0.99;
		aVel[i].m_y *= 0.99;
	}		
}

using cell = unsigned char;

void fillCells(Array<cell> &aCells, World const &world, ParticleSystem const &parSys)
{
	auto const &aPos = parSys.m_aPositions;
	auto const &aCharges = parSys.m_aCharges;
	const auto width = world.m_width;
	const auto height = world.m_height;
	Array<float> aFloatCells;
#define fcellAt(i,j) aFloatCells[((i + height) % height) * width + ((j + width) % width)]
	aFloatCells.resize(aCells.size());
	//std::cout <<aFloatCells.size() << " " << aCells.size() << "\n";
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			auto &cell = aFloatCells[i * width + j];
			cell = 0.0f;
			vec2 cellPos{float(j) + 0.5f, float(i) + 0.5f};
			for (int parId = 0; parId < aPos.size(); parId++)
			{
				auto pos = aPos[parId];
				auto charge = aCharges[parId];
   				cell += getPotential(cellPos,  pos, charge);
			}
			static int offsets[][2] =  {{1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}};
			
			for (int offsetId = 0; offsetId < 8; offsetId++)
			{
		   		for (int parId = 0; parId < aPos.size(); parId++)
				{
					auto pos = aPos[parId];
					pos.m_x += offsets[offsetId][0] * world.m_width;
					pos.m_y += offsets[offsetId][1] * world.m_height;
					cell += getPotential(cellPos, pos, aCharges[parId]);
				}
		   	}
		}
	}
	// render the boundary between positive and negative potential
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
	  	{
			cell c = ' ';
			auto val0 = fcellAt(i - 1, j);
			auto val1 = fcellAt(i + 1, j);
			auto val2 = fcellAt(i, j - 1);
			auto val3 = fcellAt(i, j + 1);
			auto self = fcellAt(i, j);
			if (self > 0.0f)
			{
				//c = '*';
			}
						
			{
				if (
					self * val0 < 0.0f && fabsf(self) < fabsf(val0)
					|| self * val1 < 0.0f && fabsf(self) < fabsf(val1)
					|| self * val2 < 0.0f && fabsf(self) < fabsf(val2)
					|| self * val3 < 0.0f && fabsf(self) < fabsf(val3)
					)
				{
					c = '#';
				}
			}
			
			aCells[i * width + j] = c;
		}
	}
	/* draw dots at the positions of the particles
	for (int parId = 0; parId < aPos.size(); parId++)
   	{
		auto pos = aPos[parId];
		int i = int(pos.m_y);
		int j = int(pos.m_x);
		i = i % height;
		j = j % width;
		aCells[i * width + j] = '@';
	}
	*/
	
}

void renderCells(Array<cell> const &aCells, World const &world)
{
	auto width = world.m_width;
	auto height = world.m_height;
	for (int i = 0; i < height; i++)
    {
		for (int j = 0; j < width; j++)
	 	{
			std::cout << aCells[i * width + j];
		}
		std::cout << "\n";
	}
}

void clear(World const &world)
{
	auto height = world.m_height;
	for (int i = 0; i <  height; i++)
	{
		std::cout << "\e[1A\r";
	}
					
}

int main(int argc, char **argv)
{
	srand(time(NULL));
	if (argc < 5 || atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atoi(argv[3]) <=0 || atof(argv[4]) <= 0.0f )
	{
		std::cout << "please provide the info <canvas width(int), canvas height(int), number of particles(int), particles charge(float)>" << std::endl;
		exit(1);
	}
	World world{atoi(argv[1]), atoi(argv[2])};
	std::cout << world.m_width << " " << world.m_height << std::endl;
	int N = atoi(argv[3]);
	ParticleSystem parSys;
	for(int i = 0; i < N; i++)
	{
		parSys.addParticle(world.genRandomPoint(), 1.0f, atof(argv[4]) * (i > N / 2 ? 1.0 : -1.0));
	}
	using cell = unsigned char;
	Array<cell> aCells;
	aCells.resize(world.m_width * world.m_height);
	renderCells(aCells, world);
	while(1)
    {
		iter(parSys, world, 1.0e-1f);
		regularize(parSys, world);
		fillCells(aCells, world, parSys);
		clear(world);
		renderCells(aCells, world);
		usleep( 10 * 1000 );
	}
    return 0;
}
