#include<iostream>
#include<particle_simulator.hpp>

struct Person{
    PS::S64    id;
    PS::S64    stat; //0:susceptible, 1:infectious, 2:recoverd/removed
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_search;
    PS::F64 t_infected;
    static inline PS::F64 time_sys;
    PS::F64vec getPos() const {
        return pos;
    }
    void setPos(const PS::F64vec & pos_new) {
        pos = pos_new;
    }
    PS::F64 getRSearch() const{
        return r_search;
    }
    void copyFromFP(const Person & fp){
	//std::cerr<<"fp.id= "<<fp.id<<" fp.stat= "<<fp.stat<<std::endl;
	id   = fp.id;
        stat = fp.stat;
        pos  = fp.pos;
	//vel  = fp.vel;
	r_search   = fp.r_search;
	//t_infected = fp.t_infected;
    }
    void copyFromForce(const Person & fo) {
	//std::cerr<<"fo.id= "<<fo.id<<" fo.stat= "<<fo.stat<<std::endl;
	//id   = fo.id;
        //pos  = fo.pos;
	//vel  = fo.vel;
	//r_search   = fo.r_search;
	stat = fo.stat;
	t_infected = fo.t_infected;
    }
    void clear() {}
    void move(const PS::F64 dt) {
	pos += vel*dt;
    }
    void dump(std::ostream & fout=std::cout){
	fout<<"id= "<<id<<" stat= "<<stat<<" pos= "<<pos<<" vel= "<<vel<<" r_search= "<<r_search<<" t_infected= "<<t_infected<<std::endl;
    }
};

void CalcInteraction(const Person pi[],
                     const PS::S32 ni,
                     const Person pj[],
                     const PS::S32 nj,
                     Person fi[]) {
    Person pj_infected[nj];
    PS::S32 nj_infected = 0;
    for(int j=0; j<nj; j++){
	if(pj[j].stat == 1){
	    pj_infected[nj_infected++] = pj[j];
	}
    }
    for(int i=0; i<ni; i++){
	if(pi[i].stat != 0){
	    fi[i].stat       = pi[i].stat;
	    fi[i].t_infected = pi[i].t_infected;
	    continue;
	}
        for(int j=0; j<nj_infected; j++){
            const auto r_crit_sq = pj_infected[j].r_search * pj_infected[j].r_search;
            const auto rij = pi[i].pos - pj_infected[j].pos;
            if(rij*rij < r_crit_sq){
                fi[i].stat = 1;
		fi[i].t_infected = Person::time_sys;
		break;
            }
        }
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
void Move(T & people, const PS::F64 dt){
    const auto n = people.getNumberOfParticleLocal();
    for(int i=0; i<n; i++){
	people[i].move(dt);
    }
}

template<typename T>
PS::S32 CountNInfected(T & people){
    const auto n = people.getNumberOfParticleLocal();
    PS::S32 ret = 0;
    for(int i=0; i<n; i++){
	if(people[i].stat==1){
	    ret++;
	}
    }
    return ret;
}

using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_SCATTER, Person, Person, Person, PS::MomentShort, PS::MomentShort, PS::SuperParticleBase, PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>;

template<typename T>
void setParticlesUniformBox(T &sys, 
                            const PS::S64 n_loc,
                            const PS::S32 seed ){
    sys.setNumberOfParticleLocal(n_loc);
    PS::F64 vol = 1.0;
    PS::F64 dis_ave = cbrt(vol/(n_loc*PS::Comm::getNumberOfProc()));
    std::cerr<<"dis_ave= "<<dis_ave<<std::endl;
    PS::MTTS mt;
    mt.init_genrand(seed);
    for(PS::S32 i = 0; i < n_loc; i++){
        sys[i].pos.x = mt.genrand_res53();
        sys[i].pos.y = mt.genrand_res53();
        sys[i].pos.z = mt.genrand_res53();
        sys[i].vel.x = 0.1*mt.genrand_res53() -0.05;
        sys[i].vel.y = 0.1*mt.genrand_res53() -0.05;
        sys[i].vel.z = 0.1*mt.genrand_res53() -0.05;
        sys[i].id    = n_loc*PS::Comm::getRank() + i;
        sys[i].r_search = dis_ave * 0.01;
    }
    Person::time_sys = 0.0;
}

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);
    PS::Initialize(argc, argv);
    PS::F64 dt = 1e-5;
    PS::F64 t_end = 10.0;
    
    PS::ParticleSystem<Person> people;
    people.initialize();
    PS::S64 n_loc = 100;
    setParticlesUniformBox(people, n_loc, PS::Comm::getRank());
    //std::cout<<"people.getNumberOfParticleLocal()= "<<people.getNumberOfParticleLocal()<<std::endl;
    //for(auto i=0; i<people.getNumberOfParticleLocal(); i++){
    //	std::cout<<"i= "<<i<<" id= "<<people[i].id<<" stat= "<<people[i].stat<<std::endl;
    //}


    
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    //dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setPosRootDomain(PS::F64vec(0.0), PS::F64vec(1.0));
    
    dinfo.decomposeDomainAll(people);
    people.adjustPositionIntoRootDomain(dinfo);
    
    people.exchangeParticle(dinfo);
    n_loc = people.getNumberOfParticleLocal();
    for(int i=0; i<n_loc; i++){
	if(people[i].id == 0){
	    people[i].stat = 1;
	    people[i].t_infected = Person::time_sys;
	}
    }
    const auto n_infected = CountNInfected(people);
    std::cout<<"time_sys= "<<Person::time_sys<<" n_infected= "<<n_infected<<std::endl;
    Tree_t tree;
    tree.initialize(n_loc);
    while(t_end > Person::time_sys){
	tree.calcForceAllAndWriteBack(CalcInteraction, people, dinfo, false);
	Move(people, dt);
	Person::time_sys += dt;
	const auto n_infected = CountNInfected(people);
	std::cout<<"time_sys= "<<Person::time_sys<<" n_infected= "<<n_infected<<std::endl;
	
	dinfo.decomposeDomainAll(people);
	people.adjustPositionIntoRootDomain(dinfo);
	people.exchangeParticle(dinfo);
    }
    PS::Finalize();
    return 0;
}
