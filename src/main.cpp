#include<iostream>
#include<particle_simulator.hpp>

struct Person
{
    PS::S64    id;
    PS::S64    stat; //0:susceptible, 1:infectious, 2:recoverd/removed
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_search;
    PS::F64 t_infected;
    static inline PS::F64 time_sys;
    static inline PS::F64vec root_length;
    PS::F64vec getPos() const 
    {
        return pos;
    }
    void setPos(const PS::F64vec & pos_new) 
    {
        pos = pos_new;
    }
    PS::F64 getRSearch() const
    {
        return r_search;
    }
    void copyFromFP(const Person & fp)
    {
	//std::cerr<<"fp.id= "<<fp.id<<" fp.stat= "<<fp.stat<<std::endl;
	id   = fp.id;
        stat = fp.stat;
        pos  = fp.pos;
	//vel  = fp.vel;
	r_search   = fp.r_search;
	t_infected = fp.t_infected;
    }
    void copyFromForce(const Person & fo) 
    {
	//std::cerr<<"fo.id= "<<fo.id<<" fo.stat= "<<fo.stat<<std::endl;
	//id   = fo.id;
        //pos  = fo.pos;
	//vel  = fo.vel;
	//r_search   = fo.r_search;
	stat = fo.stat;
	t_infected = fo.t_infected;
    }
    void clear() {}
    void move(const PS::F64 dt) 
    {
	pos += vel*dt;
    }
    void dump(std::ostream & fout=std::cout)
    {
	fout<<"id= "<<id<<" stat= "<<stat<<" pos= "<<pos<<" vel= "<<vel<<" r_search= "<<r_search<<" t_infected= "<<t_infected<<std::endl;
    }
};

void CalcInteraction(const Person pi[],
                     const PS::S32 ni,
                     const Person pj[],
                     const PS::S32 nj,
                     Person fi[]) 
{
    if(pi[0].id == 1){
	std::cerr<<"A) nj= "<<nj<<std::endl;
	for(int j=0; j<nj; j++){
	    std::cerr<<"pj[j].id= "<<pj[j].id<<" pj[j].pos= "<<pj[j].pos<<std::endl;
	}
    }

    
    for(int j=0; j<nj; j++){
	if(pj[j].id == 0){
	    std::cerr<<"B) pi[0].id= "<<pi[0].id<<" pi[0].pos= "<<pi[0].pos<<std::endl;
	}
    }
    Person pj_infected[nj];
    PS::S32 nj_infected = 0;
    for(int j=0; j<nj; j++)
    {
	if(pj[j].stat == 1)
	{
	    pj_infected[nj_infected++] = pj[j];
	}
    }
    //std::cerr<<"nj_infected= "<<nj_infected<<std::endl;
    for(int i=0; i<ni; i++)
    {
	/*
	if(pi[i].id == 34){
	    std::cerr<<"A) pi[i].stat= "<<pi[i].stat<<std::endl;
	}
	*/
	fi[i].stat       = pi[i].stat;
	fi[i].t_infected = pi[i].t_infected;	
	if(pi[i].stat != 0){
	    continue;
	}
	//std::cerr<<"pi[i].pos= "<<pi[i].pos<<std::endl;
        for(int j=0; j<nj_infected; j++)
	{
	    if(pj_infected[j].id==0){
		std::cerr<<"i= "<<i<<"pi[i].pos= "<<pi[i].pos<<std::endl;
	    }
            const auto r_crit_sq = pj_infected[j].r_search * pj_infected[j].r_search;
            const auto rij_tmp = pi[i].pos - pj_infected[j].pos;
	    auto rij = rij_tmp;
	    rij.x = (std::abs(rij.x) < (Person::root_length.x-std::abs(rij.x))) ? rij.x : ((rij.x > 0.0) ? rij.x-Person::root_length.x : rij.x+Person::root_length.x);
	    rij.y = (std::abs(rij.y) < (Person::root_length.y-std::abs(rij.y))) ? rij.y : ((rij.y > 0.0) ? rij.y-Person::root_length.y : rij.y+Person::root_length.y);
	    rij.z = (std::abs(rij.z) < (Person::root_length.z-std::abs(rij.z))) ? rij.z : ((rij.z > 0.0) ? rij.z-Person::root_length.z : rij.z+Person::root_length.z);
	    //std::cerr<<"rij= "<<rij<<std::endl;
	    //std::cerr<<"r_crit_= "<<sqrt(r_crit_sq)<<" |rij|= "<<sqrt(rij*rij)<<std::endl;
            if(rij*rij < r_crit_sq)
	    {
		//std::cerr<<"pi[i].id= "<<pi[i].id<<" pi[i].pos= "<<pi[i].pos<<" pj_infected[j].pos= "<<pj_infected[j].pos<<std::endl;
                fi[i].stat = 1;
		fi[i].t_infected = Person::time_sys;
		break;
            }
        }
	/*
	if(pi[i].id == 34){
	    std::cerr<<"B) pi[i].stat= "<<pi[i].stat<<std::endl;
	}
	*/
    }
}

/*
template<typename T, typename Func>
void applyPeople(T & people, Func func){
    const auto n = people.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
	people[i].*func();
    }
}
*/

template<typename T>
void Move(T & people, const PS::F64 dt)
{
    const auto n = people.getNumberOfParticleLocal();
    for(int i=0; i<n; i++)
    {
	people[i].move(dt);
    }
}

template<typename T>
PS::S32 CountNInfected(T & people)
{
    const auto n = people.getNumberOfParticleLocal();
    PS::S32 ret = 0;
    for(int i=0; i<n; i++)
    {
	if(people[i].stat==1)
	{
	    std::cerr<<"i= "<<i<<" people[i].id= "<<people[i].id<<std::endl;
	    ret++;
	}
    }
    return ret;
}

template<typename T>
PS::S32 CountNRemoved(T & people)
{
	const auto n = people.getNumberOfParticleLocal();
	PS::S32 ret =0;
	for(int i=0;i<n;i++)
	{
		if(people[i].stat==2)
			ret++;
	}
	return ret;
}


//using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_SCATTER, Person, Person, Person, PS::MomentShort, PS::MomentShort, PS::SuperParticleBase, PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>;
using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_SCATTER, Person, Person, Person, PS::MomentShort, PS::MomentShort, PS::SuperParticleBase>;

template<typename T>
void setParticlesUniformBox(T &sys, 
                            const PS::S64 n_loc,
                            const PS::S32 seed )
{
    sys.setNumberOfParticleLocal(n_loc);
    PS::F64 vol = 1.0;
    PS::F64 dis_ave = cbrt(vol/(n_loc*PS::Comm::getNumberOfProc()));
    std::cerr<<"dis_ave= "<<dis_ave<<std::endl;
    PS::MTTS mt;
    mt.init_genrand(seed);
    for(PS::S32 i = 0; i < n_loc; i++)
    {
        sys[i].pos.x = mt.genrand_res53();
        sys[i].pos.y = mt.genrand_res53();
	sys[i].pos.z = 0.0;
        sys[i].vel.x = 0.1*mt.genrand_res53() -0.05;
        sys[i].vel.y = 0.1*mt.genrand_res53() -0.05;
	sys[i].vel.z = 0.0;
        sys[i].id    = n_loc*PS::Comm::getRank() + i;
        sys[i].r_search = dis_ave * 0.1;
	//sys[i].r_search = 0.01;
	std::cerr<<"sys[i].r_search= "<<sys[i].r_search<<std::endl;
    }
    Person::time_sys = 0.0;
}

int main(int argc, char *argv[]) 
{
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    //PS::F64 dt = 1e-3;
    PS::F64 dt = 1e-2;
    //PS::F64 t_end = 0.1;
    PS::F64 t_end = 10.0;
   // PS::F64 time_snap = 0.0;
   // PS::S32 id_snap =0;
   FILE *gp;

   gp=popen("gnuplot -persist","w");
   fprintf(gp,"set xyplane 0\n");
   fprintf(gp,"set xrange [0.0:1.0]\n");
   fprintf(gp,"set yrange [0.0:1.0]\n");
   fprintf(gp,"set size square\n");
   fprintf(gp,"rgb(r,g,b)=65536*int(r)+256*int(g)+int(b)\n");
    
    PS::ParticleSystem<Person> people;
    people.initialize();
    PS::S64 n_loc = 1000;
    //PS::S64 n_loc = 10;
    setParticlesUniformBox(people, n_loc, PS::Comm::getRank());
    //std::cout<<"people.getNumberOfParticleLocal()= "<<people.getNumberOfParticleLocal()<<std::endl;
    //for(auto i=0; i<people.getNumberOfParticleLocal(); i++){
    //	std::cout<<"i= "<<i<<" id= "<<people[i].id<<" stat= "<<people[i].stat<<std::endl;
    //}
    
    const PS::F64 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    //dinfo.setPosRootDomain(PS::F64vec(0.0), PS::F64vec(1.0));
    dinfo.setPosRootDomainX(0.0, 1.0);
    dinfo.setPosRootDomainY(0.0, 1.0);
    dinfo.setPosRootDomainZ(0.0, 1.0);
    Person::root_length = PS::F64vec(1.0);

    dinfo.decomposeDomainAll(people);

    people.adjustPositionIntoRootDomain(dinfo);

    people.exchangeParticle(dinfo);
    n_loc = people.getNumberOfParticleLocal();
    for(int i=0; i<n_loc; i++)
    {
#if 1
	if(people[i].id < 10){
	    people[i].stat = 1;
	    people[i].t_infected = Person::time_sys;
	}
#else
	if(people[i].id == 0)
	{
	    people[i].stat = 1;
	    people[i].pos.x = 0.999;
	    people[i].pos.y = 0.0;
	    people[i].pos.z = 0.0;
	    people[i].vel.x = 0.0;
	    people[i].vel.y = 0.1;
	    people[i].vel.z = 0.0;	    
	    people[i].t_infected = Person::time_sys;
	}
	if(people[i].id == 1)
	{
	    people[i].pos.x = 0.001;
	    people[i].pos.y = 0.05;
	    people[i].pos.z = 0.0;
	    people[i].vel.x = 0.0;
	    people[i].vel.y = -0.1;
	    people[i].vel.z = 0.0;	    
	    people[i].t_infected = Person::time_sys;
	}
#endif
    }
    const auto n_infected = CountNInfected(people);
    const auto n_removed = CountNRemoved(people);

    std::cout<<"|time_sys = "<<Person::time_sys<<" |n_infected = "<<n_infected<<std::endl;
    Tree_t tree;
    //tree.initialize(n_loc, 0.5, 1, 1);
    tree.initialize(n_loc);
    while(t_end > Person::time_sys)
    {


	tree.calcForceAllAndWriteBack(CalcInteraction, people, dinfo, false);
	/*
	for(int i=0; i<people.getNumberOfParticleLocal(); i++)	{
	    if(people[i].id == 34){
		std::cerr<<"C) people[i].stat= "<<people[i].stat<<std::endl;
	    }
	}
	*/
	Move(people, dt);
	Person::time_sys += dt;
	const auto n_infected = CountNInfected(people);
	const auto n_removed = CountNRemoved(people);
	
	for(int i=0;i<n_loc; i++)
	{
	    if(people[i].stat == 1 && (Person::time_sys - people[i].t_infected) > 3.0)
		//if(0)
		{
			people[i].stat = 2;

		}
	}
	std::cout<<"|time_sys = "<<Person::time_sys<<" |n_infected = "<<n_infected<<" |n_removed = "<<n_removed<<std::endl;
	
	dinfo.decomposeDomainAll(people);

	people.adjustPositionIntoRootDomain(dinfo);

	people.exchangeParticle(dinfo);

      //fprintf(gp, "plot '-' using 1:2:3:(rgb($4,$5,$6)) w p pt 7 ps 0.8 lc rgb variable\n");
      fprintf(gp, "plot '-' using 1:2:(rgb($4,$5,$6)) w p pt 7 ps 0.8 lc rgb variable\n");
    	for(PS::S32 i = 0; i < n_loc; i++)
	{
		fprintf(gp,"%f, %f, %f ",people[i].pos.x, people[i].pos.y, people[i].pos.z);
		if(people[i].stat == 1)
			fprintf(gp,"255 0 0 \n");
		else if(people[i].stat == 2)
			fprintf(gp,"0 0 255 \n");
		else
			fprintf(gp,"0 255 0 \n");

	}
	fprintf(gp, "e\n");
    }
    pclose(gp);
    PS::Finalize();
    return 0;
}
