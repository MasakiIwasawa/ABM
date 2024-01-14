#include <cstdlib>
#include <iostream>
#include <particle_simulator.hpp>
#include <random>
#define DEBUG_LEVEL 0
#define INDIVISUAL_WALL
#define DENSITY_CONTRAST
constexpr int SEED = 0;
std::mt19937 RND_ENG(SEED);

constexpr PS::F64 DOMAIN_LENGTH = 1.0;
constexpr PS::F64 DOMAIN_AREA   = DOMAIN_LENGTH * DOMAIN_LENGTH;
constexpr PS::S64 N_PEOPLE_GLB = 1000;
constexpr PS::S32 n_infected_init = 3;
constexpr PS::F64 DENS_IN_DENS_REGION_NOR = 100.0;
//constexpr PS::F64 DT_INFECTED_NOR = 3.0; // infected duration normalized by t_coll
constexpr PS::F64 DT_INFECTED_NOR = 0.3; // infected duration normalized by t_coll
constexpr PS::F64 DT_RECOVERD_NOR = 1.0; // recoverd duration normalized by t_coll
constexpr PS::F64 LENGTH_DENS_REGION = 0.01;
constexpr PS::S32 N_DENS_REGION_1D = 2;
constexpr PS::F64 VEL_AVE = 0.05;
constexpr PS::F64 R_INFECTED_NOR = 0.05; // infected radius normalized by average particle distance

constexpr auto DENS_GLB = (PS::F64)N_PEOPLE_GLB / DOMAIN_AREA;
constexpr auto R_INFECTED = R_INFECTED_NOR * sqrt(DOMAIN_AREA / N_PEOPLE_GLB);
constexpr auto SIGMA = R_INFECTED * 2.0; // cross section

void Reflect(PS::F64 &x, PS::F64 &v, PS::F64 x_min, PS::F64 x_max) {
    auto x_tmp = std::clamp(x, x_min, x_max);
    v = (x_tmp == x) ? v : -v;
    x = x_tmp + (x_tmp - x);
}
void Reflect(PS::F64vec &x, PS::F64vec &v, const PS::F64vec x_min,
             const PS::F64vec x_max) {
    Reflect(x.x, v.x, x_min.x, x_max.x);
    Reflect(x.y, v.y, x_min.y, x_max.y);
    Reflect(x.z, v.z, x_min.z, x_max.z);
}
PS::F64vec GetRandomPos(const PS::F64ort &bnd) {
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    const auto x = rand(RND_ENG);
    const auto y = rand(RND_ENG);
    const auto x_new = bnd.low_.x + (bnd.high_.x - bnd.low_.x) * x;
    const auto y_new = bnd.low_.y + (bnd.high_.y - bnd.low_.y) * y;
    const auto z_new = 0.0;
    return PS::F64vec(x_new, y_new, z_new);
}

// 人口密度を変更する.
template <typename Tsys> void ChangeDensity(Tsys &people) {
    //const PS::F64 dens_glb = N_PEOPLE_GLB / DOMAIN_AREA;
    PS::S32 n_dens_region_2d = N_DENS_REGION_1D * N_DENS_REGION_1D;
    PS::F64ort pos_dens[N_DENS_REGION_1D][N_DENS_REGION_1D];
    constexpr PS::F64 dx_half = LENGTH_DENS_REGION;
    const PS::F64 area_dens_region = (dx_half * 2.0) * (dx_half * 2.0);
    const PS::F64 dens_in_dens_resion = DENS_GLB * DENS_IN_DENS_REGION_NOR;
    const PS::S32 n_in_dens_region = dens_in_dens_resion * area_dens_region * n_dens_region_2d;
    std::cerr<<"n_in_dens_region= "<<n_in_dens_region<<std::endl;
    for (int i = 0; i < N_DENS_REGION_1D; i++) {
        for (int j = 0; j < N_DENS_REGION_1D; j++) {
            PS::F64 x_cen = (1.0 / (N_DENS_REGION_1D + 1)) * (i + 1);
            PS::F64 y_cen = (1.0 / (N_DENS_REGION_1D + 1)) * (j + 1);
            pos_dens[i][j].low_ =
                PS::F64vec(x_cen - dx_half, y_cen - dx_half, 0.0);
            pos_dens[i][j].high_ =
                PS::F64vec(x_cen + dx_half, y_cen + dx_half, 0.0);
            std::cerr << "i= " << i << " j= " << j
                      << " pos_dens[i][j]= " << pos_dens[i][j] << std::endl;
        }
    }
    const auto avg_particle_distance = sqrt(1.0 / DENS_GLB);
    const auto half_wall_length = avg_particle_distance * 2.0;
    const PS::F64vec half_wall_length_vec(half_wall_length);
    const auto n_loc = people.getNumberOfParticleLocal();
    std::uniform_int_distribution<> rand(0, N_DENS_REGION_1D - 1);
    for (int i = 0; i < n_in_dens_region; i++) {
        const auto i0 = rand(RND_ENG);
        const auto j0 = rand(RND_ENG);
        const auto pos0 = GetRandomPos(pos_dens[i0][j0]);
        people[i].pos = pos0;
        const auto pos_low =
            (people[i].pos - half_wall_length_vec).applyEach([](const auto x) {
                return std::max(x, 0.0);
            });
        const auto pos_high =
            (people[i].pos + half_wall_length_vec).applyEach([](const auto x) {
                return std::min(x, 1.0);
            });
        people[i].wall = PS::F64ort(pos_low, pos_high);
    }
}

struct Person {
    PS::S64 id;
    PS::S64 stat; // 0:susceptible, 1:infectious, 2:recoverd/removed
    PS::F64vec pos;
    PS::F64vec vel;
    PS::F64 r_search;
    PS::F64 t_infected;
    int group_id;
    PS::F64ort wall;
    static inline PS::F64 time_sys;
    static inline PS::F64vec root_length;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec &pos_new) { pos = pos_new; }
    PS::F64 getRSearch() const { return r_search; }
    void copyFromFP(const Person &fp) {
        // std::cerr<<"fp.id= "<<fp.id<<" fp.stat= "<<fp.stat<<std::endl;
        id = fp.id;
        stat = fp.stat;
        pos = fp.pos;
        // vel  = fp.vel;
        r_search = fp.r_search;
        t_infected = fp.t_infected;
    }
    void copyFromForce(const Person &fo) {
        // std::cerr<<"fo.id= "<<fo.id<<" fo.stat= "<<fo.stat<<std::endl;
        // id   = fo.id;
        // pos  = fo.pos;
        // vel  = fo.vel;
        // r_search   = fo.r_search;
        stat = fo.stat;
        t_infected = fo.t_infected;
    }
    void clear() {}
    void move(const PS::F64 dt) { pos += vel * dt; }
    void update(const PS::F64 dt, const int mode = 0) {
        if (mode == 0) {
            move(dt);
        } else if (mode == 1) {
            pos += vel * dt;
            Reflect(pos, vel, wall.low_, wall.high_);
        } else {
            std::cerr << "mode= " << mode << " is not supported." << std::endl;
            PS::Abort();
        }
    }
    void set_wall(const PS::F64ort & w = PS::F64ort(PS::F64vec(0.0), PS::F64vec(1.0, 1.0, 0.0))){
        wall.low_.x  = w.low_.x;
        wall.low_.y  = w.low_.y;
        wall.low_.z  = w.low_.z;
        wall.high_.x = w.high_.x;
        wall.high_.y = w.high_.y;
        wall.high_.z = w.high_.z;
    }
    void dump(std::ostream &fout = std::cout) {
        fout << "id= " << id << " stat= " << stat << " pos= " << pos
             << " vel= " << vel << " r_search= " << r_search
             << " t_infected= " << t_infected << std::endl;
    }
};

void CalcInteraction(const Person pi[], const PS::S32 ni, const Person pj[],
                     const PS::S32 nj, Person fi[]) {
    /*
    if(pi[0].id == 1){
    std::cerr<<"A) nj= "<<nj<<std::endl;
    for(int j=0; j<nj; j++){
        std::cerr<<"pj[j].id= "<<pj[j].id<<" pj[j].pos= "<<pj[j].pos<<std::endl;
    }
    }
    */

    /*
    for(int j=0; j<nj; j++){
    if(pj[j].id == 0){
        std::cerr<<"B) pi[0].id= "<<pi[0].id<<" pi[0].pos=
    "<<pi[0].pos<<std::endl;
    }
    }
    */
    Person pj_infected[nj];
    PS::S32 nj_infected = 0;
    for (int j = 0; j < nj; j++) {
        if (pj[j].stat == 1) {
            pj_infected[nj_infected++] = pj[j];
        }
    }
    // std::cerr<<"nj_infected= "<<nj_infected<<std::endl;
    for (int i = 0; i < ni; i++) {
        /*
        if(pi[i].id == 34){
            std::cerr<<"A) pi[i].stat= "<<pi[i].stat<<std::endl;
        }
        */
        fi[i].stat = pi[i].stat;
        fi[i].t_infected = pi[i].t_infected;
        if (pi[i].stat != 0) {
            continue;
        }
        // std::cerr<<"pi[i].pos= "<<pi[i].pos<<std::endl;
        for (int j = 0; j < nj_infected; j++) {
            // if(pj_infected[j].id==0){
            //	std::cerr<<"i= "<<i<<"pi[i].pos= "<<pi[i].pos<<std::endl;
            // }
            const auto r_crit_sq =
                pj_infected[j].r_search * pj_infected[j].r_search;
            const auto rij_tmp = pi[i].pos - pj_infected[j].pos;
            auto rij = rij_tmp;
            rij.x =
                (std::abs(rij.x) < (Person::root_length.x - std::abs(rij.x)))
                    ? rij.x
                    : ((rij.x > 0.0) ? rij.x - Person::root_length.x
                                     : rij.x + Person::root_length.x);
            rij.y =
                (std::abs(rij.y) < (Person::root_length.y - std::abs(rij.y)))
                    ? rij.y
                    : ((rij.y > 0.0) ? rij.y - Person::root_length.y
                                     : rij.y + Person::root_length.y);
            rij.z =
                (std::abs(rij.z) < (Person::root_length.z - std::abs(rij.z)))
                    ? rij.z
                    : ((rij.z > 0.0) ? rij.z - Person::root_length.z
                                     : rij.z + Person::root_length.z);
            // std::cerr<<"rij= "<<rij<<std::endl;
            // std::cerr<<"r_crit_= "<<sqrt(r_crit_sq)<<" |rij|=
            // "<<sqrt(rij*rij)<<std::endl;
            if (rij * rij < r_crit_sq) {
                // std::cerr<<"pi[i].id= "<<pi[i].id<<" pi[i].pos=
                // "<<pi[i].pos<<" pj_infected[j].pos=
                // "<<pj_infected[j].pos<<std::endl;
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

template <typename T> void Move(T &people, const PS::F64 dt) {
    const auto n = people.getNumberOfParticleLocal();
    for (int i = 0; i < n; i++) {
#if defined(INDIVISUAL_WALL)
        people[i].update(dt, 1);
#else
        people[i].update(dt, 0);
#endif
        // people[i].move(dt);
    }
}

template <typename T> PS::S32 CountNInfected(T &people) {
    const auto n = people.getNumberOfParticleLocal();
    PS::S32 ret = 0;
    for (int i = 0; i < n; i++) {
        if (people[i].stat == 1) {
            // std::cerr<<"i= "<<i<<" people[i].id= "<<people[i].id<<std::endl;
            ret++;
        }
    }
    return ret;
}

template <typename T> PS::S32 CountNRemoved(T &people) {
    const auto n = people.getNumberOfParticleLocal();
    PS::S32 ret = 0;
    for (int i = 0; i < n; i++) {
        if (people[i].stat == 2)
            ret++;
    }
    return ret;
}

// using Tree_t = PS::TreeForForce<PS::SEARCH_MODE_SCATTER, Person, Person,
// Person, PS::MomentShort, PS::MomentShort, PS::SuperParticleBase,
// PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>;
using Tree_t =
    PS::TreeForForce<PS::SEARCH_MODE_SCATTER, Person, Person, Person,
                     PS::MomentShort, PS::MomentShort, PS::SuperParticleBase>;

template <typename T>
void setParticlesUniformBox(T &sys, const PS::S64 n_loc, const PS::S32 seed = 0) {
    sys.setNumberOfParticleLocal(n_loc);
    PS::F64 area = 1.0 * 1.0;
    PS::F64 dens = (n_loc * PS::Comm::getNumberOfProc()) / DOMAIN_AREA;
    PS::F64 dis_ave = sqrt(DOMAIN_AREA / (n_loc * PS::Comm::getNumberOfProc()));
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    for (PS::S32 i = 0; i < n_loc; i++) {
        sys[i].pos.x = rand(RND_ENG);
        sys[i].pos.y = rand(RND_ENG);
        sys[i].pos.z = 0.0;
        sys[i].vel.x = 0.1 * rand(RND_ENG) - 0.05;
        sys[i].vel.y = 0.1 * rand(RND_ENG) - 0.05;
        sys[i].vel.z = 0.0;
        sys[i].id = n_loc * PS::Comm::getRank() + i;
        sys[i].r_search = dis_ave * R_INFECTED_NOR;
        sys[i].set_wall();
        // std::cerr<<"sys[i].r_search= "<<sys[i].r_search<<std::endl;
    }
    // const auto vel_ave = 0.05;
    // const auto Tcoll = 1.0/(dens * sys[0].r_search * vel_ave);
    // std::cerr<<"Tcoll= "<<Tcoll<<std::endl;
    // const auto Tcross = sys[0].r_search / vel_ave;
    // std::cerr<<"Tcross= "<<Tcross<<std::endl;
    // exit(1);
    Person::time_sys = 0.0;
}

template <typename T>
void wall(T &people, const PS::S64 n_loc, const PS::F64 dt) {
    // const auto wall1_l,wall1_r,wall1_o,wall1_u;
    // const auto wall2_l,wall2_r,wall2_o,wall2_u;
    const auto wall = 0.5;
    Person hist[n_loc];
    for (int i = 0; i < n_loc; i++) {
        hist[i] = people[i];
        hist[i].pos = people[i].pos - people[i].vel * dt;
        if (hist[i].pos.y > wall && people[i].pos.y <= wall) {
            people[i].vel.y = -people[i].vel.y;
            people[i].pos.y = wall + (wall - people[i].pos.y);
#if DEBUG_LEVEL >= 1
            if (hist[i].pos.y <= wall || people[i].pos.y <= wall) {
                std::cout << "error1" << std::endl;
                std::cout << "i= " << i << " pos= " << people[i].pos
                          << " vel=" << people[i].vel << std::endl;
                std::cout << "i= " << i << " pos= " << hist[i].pos
                          << " vel=" << hist[i].vel << std::endl;
            }
#endif
        } else if (hist[i].pos.y < wall && people[i].pos.y >= wall) {
            people[i].vel.y = -people[i].vel.y;
            people[i].pos.y = wall + (wall - people[i].pos.y);
#if DEBUG_LEVEL >= 1
            if (hist[i].pos.y >= wall || people[i].pos.y >= wall) {
                std::cout << "error2" << std::endl;
                std::cout << "i= " << i << " pos= " << people[i].pos
                          << " vel=" << people[i].vel << std::endl;
                std::cout << "i= " << i << " pos= " << hist[i].pos
                          << " vel=" << hist[i].vel << std::endl;
            }
#endif
        }
    }
    // for debug
#if DEBUG_LEVEL >= 1
    int n_up_loc = 0;
    int n_dw_loc = 0;
    for (int i = 0; i < n_loc; i++) {
        if (people[i].pos.y >= wall) {
            n_up_loc++;
        } else {
            n_dw_loc++;
        }
    }
    int n_up = PS::Comm::getSum(n_up_loc);
    int n_dw = PS::Comm::getSum(n_dw_loc);
    std::cout << "n_up= " << n_up << " n_dw= " << n_dw << std::endl;
#endif
}

int main(int argc, char *argv[]) {
    std::cout << std::setprecision(15);
    std::cerr << std::setprecision(15);
    PS::Initialize(argc, argv);
    const auto my_rank = PS::Comm::getRank();
    const auto n_proc  = PS::Comm::getNumberOfProc();
    FILE *gp;
    gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set xyplane 0\n");
    fprintf(gp, "set xrange [0.0:1.0]\n");
    fprintf(gp, "set yrange [0.0:1.0]\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "rgb(r,g,b)=65536*int(r)+256*int(g)+int(b)\n");

    PS::ParticleSystem<Person> people;
    people.initialize();
    PS::S64 n_loc = N_PEOPLE_GLB / n_proc;
    n_loc = (my_rank == n_proc-1) ? N_PEOPLE_GLB - n_loc * (n_proc - 1) : n_loc;
    setParticlesUniformBox(people, n_loc, PS::Comm::getRank());
#if defined(DENSITY_CONTRAST)
    ChangeDensity(people);
#endif
    const auto n_glb = people.getNumberOfParticleGlobal();
    //const auto area = 1.0;
    //const auto vel_ave = 0.05;
    //const auto dens = (PS::F64)n_glb / area;
    //const auto sigma = people[0].r_search * 2.0;
    const auto t_coll = 1.0 / (DENS_GLB * SIGMA * VEL_AVE);
    const auto dt_infected  = DT_INFECTED_NOR * t_coll;
    const auto dt_recoverd = DT_RECOVERD_NOR * t_coll;    
    // const auto t_end = t_coll * 10.0;
    const auto t_end = t_coll * 300.0;
    //const auto t_recover = t_coll * 3.0;
    const auto dt = SIGMA / (VEL_AVE * 2.0) * 0.1;
    std::cout << "# n_glb= " << n_glb << " VEL_AVE= " << VEL_AVE
              << " DENS_GLB= " << DENS_GLB << " SIGMA= " << SIGMA
              << " t_coll= " << t_coll << " t_end= " << t_end
              << " dt_infected= " << dt_infected
              << " dt_recoverd= " << dt_recoverd
              << " dt= " << dt << std::endl;
    n_loc = people.getNumberOfParticleLocal();
    
    // exit(1);
    // std::cout<<"people.getNumberOfParticleLocal()=
    // "<<people.getNumberOfParticleLocal()<<std::endl; for(auto i=0;
    // i<people.getNumberOfParticleLocal(); i++){
    //	std::cout<<"i= "<<i<<" id= "<<people[i].id<<" stat=
    //"<<people[i].stat<<std::endl;
    // }
    // int xx[n_loc]={0},yy[n_loc]={0},zz[n_loc]={0};
    // int xx,yy,zz;

    const PS::F64 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    // dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    // dinfo.setPosRootDomain(PS::F64vec(0.0), PS::F64vec(1.0));
    dinfo.setPosRootDomainX(0.0, 1.0);
    dinfo.setPosRootDomainY(0.0, 1.0);
    dinfo.setPosRootDomainZ(0.0, 1.0);
    Person::root_length = PS::F64vec(1.0);

    dinfo.decomposeDomainAll(people);

    people.adjustPositionIntoRootDomain(dinfo);

    people.exchangeParticle(dinfo);
    n_loc = people.getNumberOfParticleLocal();
    for (int i = 0; i < n_loc; i++) {
#if 1
        if (people[i].id < n_infected_init) {
            people[i].stat = 1;
            people[i].t_infected = Person::time_sys;
        }
#else
        if (people[i].id == 0) {
            people[i].stat = 1;
            people[i].pos.x = 0.999;
            people[i].pos.y = 0.0;
            people[i].pos.z = 0.0;
            people[i].vel.x = 0.0;
            people[i].vel.y = 0.1;
            people[i].vel.z = 0.0;
            people[i].t_infected = Person::time_sys;
        }
        if (people[i].id == 1) {
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

    std::cout << "|time_sys = " << Person::time_sys
              << " |n_infected = " << n_infected << std::endl;
    Tree_t tree;
    // tree.initialize(n_loc, 0.5, 1, 1);
    tree.initialize(n_loc);
    PS::S64 n_loop = 0;
    while (t_end > Person::time_sys) {
        /*
            xx=0;
            yy=0;
            zz=0;
        */

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

        for (int i = 0; i < n_loc; i++) {
            if (people[i].stat == 1 &&
                (Person::time_sys - people[i].t_infected) > dt_infected)
            // if(0)
            {
                people[i].stat = 2;
            } else if (people[i].stat == 2 &&
                       (Person::time_sys - (people[i].t_infected + dt_infected)) >
                           dt_recoverd) {
                people[i].stat = 0;
            }
        }
        // std::cout<<" r_search= "<<people[0].r_search<<std::endl;
        std::cout << "|time_sys = " << Person::time_sys
                  << " |n_infected = " << n_infected
                  << " |n_removed = " << n_removed << " |other = "
                  << (int)n_loc - (int)n_infected - (int)n_removed << " |"
                  << std::endl;
        // fprintf(gp, "set term qt 2\n");
        // fprintf(gp,"%d,%d,%d \n",(int)n_loc-(int)n_infected-(int)n_removed,
        // n_infected, n_removed);

        dinfo.decomposeDomainAll(people);

        // wall(people,n_loc,dt);

        people.adjustPositionIntoRootDomain(dinfo);

        people.exchangeParticle(dinfo);

        if (n_loop % 100 == 0) {
            // fprintf(gp, "plot '-' using 1:2:3:(rgb($4,$5,$6)) w p pt 7 ps 0.8
            // lc rgb variable\n"); fprintf(gp, "set term qt 1\n");
            fprintf(gp, "plot '-' using 1:2:(rgb($4,$5,$6)) w p pt 7 ps 1.0 lc "
                        "rgb variable\n");
            for (PS::S32 i = 0; i < n_loc; i++) {
                fprintf(gp, "%f, %f, %f ", people[i].pos.x, people[i].pos.y,
                        people[i].pos.z);
                if (people[i].stat == 1) {
                    fprintf(gp, "255 0 0 \n");
                    // xx[i]++;
                } else if (people[i].stat == 2) {
                    fprintf(gp, "0 0 255 \n");
                    // yy[i]++;
                } else {
                    fprintf(gp, "0 255 0 \n");
                    // zz[i]++;
                }
                // fprintf(gp,"set term qt 2\n");
                // fprintf(gp,"plot (int)n_loc-(int)n_infected-(int)n_removed,
                // n_infected, n_removed\n");
            }
            fprintf(gp, "e\n");
            // exit(1);
        }
        n_loop++;
    }

    pclose(gp);
    PS::Finalize();
    return 0;
}
