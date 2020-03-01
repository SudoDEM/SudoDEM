#pragma once

#include<sudodem/core/TimeStepper.hpp>

class Integrator;


typedef std::vector<Real> stateVector;// Currently, we are unable to use Eigen library within odeint

/*Observer used to update the state of the scene*/
class observer
{
	Integrator* integrator;
public:
	observer(Integrator* _in):integrator(_in){}
        void operator()( const stateVector& /* x */ , Real /* t */ ) const;
};

//[ ode_wrapper
template< class Obj , class Mem >
class ode_wrapper
{
    Obj m_obj;
    Mem m_mem;

public:

    ode_wrapper( Obj obj , Mem mem ) : m_obj( obj ) , m_mem( mem ) { }

    template< class State , class Deriv , class Time >
    void operator()( const State &x , Deriv &dxdt , Time t )
    {
        (m_obj.*m_mem)( x , dxdt , t );
    }
};

template< class Obj , class Mem >
ode_wrapper< Obj , Mem > make_ode_wrapper( Obj obj , Mem mem )
{
    return ode_wrapper< Obj , Mem >( obj , mem );
}
//]



class Integrator: public TimeStepper {

		public:

		stateVector accumstateofthescene;//pos+vel

		stateVector accumstatedotofthescene;//only the accelerations

		stateVector resetstate;//last state before integration attempt

		Real timeresetvalue;

		inline void evaluateQuaternions(const stateVector &); //evaluate quaternions after integration

		typedef vector<vector<shared_ptr<Engine> > > slaveContainer;

		#ifdef SUDODEM_OPENMP
			vector<Real> threadMaxVelocitySq;
		#endif

		virtual void action();

		virtual void system(const stateVector&, stateVector&, Real); //System function to calculate the derivatives of states

		virtual bool isActivated(){return true;}
		// py access
		boost::python::list slaves_get();

		stateVector& getSceneStateDot();

		bool saveCurrentState(Scene const* ourscene);//Before any integration attempt state of the scene should be saved.

		bool resetLastState(void);//Before any integration attempt state of the scene should be saved.

		void slaves_set(const boost::python::list& slaves);

		stateVector& getCurrentStates(void);

		bool setCurrentStates(stateVector);

		Real updatingDispFactor;//(experimental) Displacement factor used to trigger bound update: the bound is updated only if updatingDispFactor*disp>sweepDist when >0, else all bounds are updated.

		void saveMaximaDisplacement(const shared_ptr<Body>& b);

		#ifdef SUDODEM_OPENMP
			void ensureSync(); bool syncEnsured;
		#endif


		SUDODEM_CLASS_BASE_DOC_ATTRS_CTOR_PY(Integrator,TimeStepper,"Integration Engine Interface.",
		((slaveContainer,slaves,,,"[will be overridden]"))
		((Real,integrationsteps,,,"all integrationsteps count as all succesfull substeps"))
		((Real,maxVelocitySq,NaN,,"store square of max. velocity, for informative purposes; computed again at every step. |yupdate|"))
		,
		/*ctor*/
		#ifdef SUDODEM_OPENMP
			threadMaxVelocitySq.resize(omp_get_max_threads()); syncEnsured=false;
		#endif
		,
		/*py*/

		.add_property("slaves",&Integrator::slaves_get,&Integrator::slaves_set,"List of lists of Engines to calculate the force acting on the particles;  to obtain the derivatives of the states, engines inside will be run sequentially.");
	);
};
REGISTER_SERIALIZABLE(Integrator);


