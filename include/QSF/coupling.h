/* ::::::::::::::::::::::::::::::: Gauge types :::::::::::::::::::::::::::::::::
Feel free to add more gauges recognized by the physics community. 			  */
struct VelocityGauge { //p.A
	static constexpr std::string_view CouplingName = "VG";
};
struct LengthGauge { // r.F
	static constexpr std::string_view CouplingName = "LG";
};
struct CoulombGauge { // Full coupling: p^2/2 - p*A + A^2/2
	static constexpr std::string_view CouplingName = "CG";
};
struct MixedGauge { //p.A+r.F
	static constexpr std::string_view CouplingName = "MG";
};

/* ::::::::::::::::::::::::::::::: Approx model ::::::::::::::::::::::::::::::::
Beside gauge type one can use dipole approximation or reduced dimensionality
models such as models, eg. EckhardtSachaModel which assumes one electron per
axis in 2D/3D and DipoleApprox(imation). */
struct NoApprox {
	static constexpr std::string_view CouplingName = "";
};
struct DipoleApprox {
	static constexpr std::string_view CouplingName = "DA";
};
struct EckhardSachaModel : DipoleApprox {
	static constexpr std::string_view CouplingName = "ES";
};

/* :::::::::::::::::::::::::::::::: Couplings ::::::::::::::::::::::::::::::::::
Important: Gauge types are generic types which do not include any coordinate x/p
multipliers, as gauges act diffently for different Hamiltonians (if exist).
Hence the multipliers are applied conditionally at the Hamiltonian level.
Coupling type does two things:
1) gathers fields acting on different axesand converts them to potentials
as needed */
template <class ApproxModel, class ... Fields>
struct CouplingBase;

// template <>
// struct CouplingBase<DipoleApprox> : _XBUFFER
// {
// 	// static_assert(is_unique<Fields...>, "All Fields need to act on different axes!");

// 	// tuple<Fields...> fields;
// 	// double last_values[sizeof...(Fields)]{ 0 };
// 	// double last_values_backup[sizeof...(Fields)]{ 0 };
// 	CouplingBase(Section& settings) {}
// 	static double maxPulseDuration() { return 0; }
// 	// void reset() { last_values{ 0 }; last_values_backup{ 0 }; }
// 	double operator[](AXIS ax)
// 	{
// 		return 0;
// 	}
// };


template <class ... Fields>
struct CouplingBase<DipoleApprox, Fields...>
{
	static_assert(is_unique<Fields...>, "All Fields need to act on different axes!");
	static constexpr uind size = sizeof...(Fields);
	std::tuple<Fields...> fields;
	double last_values[Max(sizeof...(Fields), 1)]{ 0 };
	double last_values_backup[Max(sizeof...(Fields), 1)]{ 0 };

	CouplingBase(Section& settings) : fields(Fields{ settings }...) {}
	CouplingBase(Fields...fields) :fields(fields...) {}
	static double maxPulseDuration() { return Max(Fields::maxPulseDuration()...); }

	// void reset() { last_values{ 0 }; last_values_backup{ 0 }; }
	double operator[](AXIS ax)
	{
		return last_values[int(ax)];
	}
};

template <class ... Fields>
struct CouplingBase<NoApprox, Fields...>
{
	std::tuple<Fields...> fields;
	double middle_values[sizeof...(Fields)]{ 0 };
	double middle_values_backup[sizeof...(Fields)]{ 0 };

	static double maxPulseDuration() { return Max(Fields::maxPulseDuration()...); }
	// void reset() { middle_values{ 0 }; middle_values_backup{ 0 }; }

	double operator[](AXIS ax)
	{
		return middle_values[ax];
	}
};

template <class GaugeType, class ... Fields>
struct DipoleCoupling;

template <class ... Fields>
struct DipoleCoupling<VelocityGauge, Fields...> : CouplingBase<DipoleApprox, Fields...>
{
	using base = CouplingBase<DipoleApprox, Fields...>;
	static constexpr REP couplesInRep = REP::P;

	double lastTime;
	DipoleCoupling(Section& settings) : base(settings) {}

	DipoleCoupling(Fields...fields) :base(fields...) {}

	DipoleCoupling(const DipoleCoupling&) = default;
	// A(t) = -int_0^t E(t) dt
	template <REP R, OPTIMS opt>
	inline void precalc() const
	{
		if constexpr (R == REP::P)
		{
			// int i = 0;
			// ([&]
			//  {
			// 	 base::last_values[i++] -= (timer - lastTime) * (std::get <Fields>(base::fields).operator()(timer));
			//  }(), ...);
			// lastTime = timer;
		}
	}
	inline double operator()()
	{

	}
	// using base::fields;
};


// template<typename ... Args>
// struct CouplingOperator : COMPUTATION<double, Args...>, _XBUFFER {
// 	using XField = NoField;
// 	using YField = NoField;
// 	using ZField = NoField;
// 	static void backup() {	}
// 	static void restore() {	}
// 	inline static constexpr double lastXValue = 0;
// 	inline static constexpr double lastYValue = 0;
// 	inline static constexpr double lastZValue = 0;
// 	static double maxPulseDuration() { return 0.0; }
// };

// template<typename FieldX, typename FieldY, typename FieldZ>
// struct CouplingOperator<FieldX, FieldY, FieldZ> : COMPUTATION<double, FieldX, FieldY, FieldZ>, _XBUFFER
// {
// private:
// 	inline static double lastXValue_copy = 0;
// 	inline static double lastYValue_copy = 0;
// 	inline static double lastZValue_copy = 0;
// public:
// 	using XField = FieldX;
// 	using YField = FieldY;
// 	using ZField = FieldZ;
// 	inline static double lastXValue = 0;
// 	inline static double lastYValue = 0;
// 	inline static double lastZValue = 0;
// 	static void backup()
// 	{
// 		lastXValue_copy = lastXValue;
// 		lastYValue_copy = lastYValue;
// 		lastZValue_copy = lastZValue;

// 	}
// 	static void restore()
// 	{
// 		lastXValue = lastXValue_copy;
// 		lastYValue = lastYValue_copy;
// 		lastZValue = lastZValue_copy;
// 	}
// 	static double maxPulseDuration()
// 	{
// 		double x_dur = 0, y_dur = 0, z_dur = 0;
// 		x_dur = XField::instance.maxPulseDuration();
// 		y_dur = YField::instance.maxPulseDuration();
// 		z_dur = ZField::instance.maxPulseDuration();
// 		return (x_dur > y_dur ? x_dur : (y_dur > z_dur ? y_dur : z_dur));
// 	}
// 	//CALC
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		lastXValue = 0.0;
// 		lastYValue = 0.0;
// 		lastZValue = 0.0;
// 	}
// 	template <REP R, OPTIMS opt, typename F>
// 	double calc() const
// 	{
// 		if constexpr (is_same_v<F, XField>) return lastXValue;
// 		if constexpr (is_same_v<F, YField>) return lastYValue;
// 		else return lastZValue;

// 	}

// 	template <REP R, typename WHEN> static constexpr inline bool canRun =
// 		COMPUTATION<double, FieldX, FieldY, FieldZ>::template goodRep<R> && is_same_v<WHEN, LATE>;
// };


// template<typename FieldX, typename FieldY>
// struct CouplingOperator<FieldX, FieldY> : COMPUTATION<double, FieldX, FieldY>, _XBUFFER
// {
// private:
// 	inline static double lastXValue_copy = 0;
// 	inline static double lastYValue_copy = 0;
// public:
// 	using XField = FieldX;
// 	using YField = FieldY;
// 	using ZField = NoField;
// 	inline static double lastXValue = 0;
// 	inline static double lastYValue = 0;
// 	inline static constexpr double lastZValue = 0;
// 	static void backup()
// 	{
// 		lastXValue_copy = lastXValue;
// 		lastYValue_copy = lastYValue;

// 	}
// 	static void restore()
// 	{
// 		lastXValue = lastXValue_copy;
// 		lastYValue = lastYValue_copy;
// 	}
// 	static double maxPulseDuration()
// 	{
// 		double x_dur = 0, y_dur = 0;
// 		x_dur = XField::instance.maxPulseDuration();
// 		y_dur = YField::instance.maxPulseDuration();

// 		return (x_dur > y_dur ? x_dur : y_dur);
// 	}
// 	//CALC
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		lastXValue = 0.0;
// 		lastYValue = 0.0;
// 	}
// 	template <REP R, OPTIMS opt, typename F>
// 	double calc() const
// 	{
// 		if constexpr (is_same_v<F, XField>) return lastXValue;
// 		else return lastYValue;
// 	}
// 	template <REP R, typename WHEN> static constexpr inline bool canRun =
// 		COMPUTATION<double, FieldX, FieldY>::template goodRep<R> && is_same_v<WHEN, LATE>;
// };
// template<typename FieldX>
// struct CouplingOperator<FieldX> : COMPUTATION<double, FieldX>, _XBUFFER
// {
// private:
// 	inline static double lastXValue_copy = 0;
// 	inline static double lastYValue_copy = 0;
// 	inline static double lastZValue_copy = 0;
// public:
// 	using XField = FieldX;
// 	using YField = NoField;
// 	using ZField = NoField;
// 	inline static double lastXValue = 0;
// 	inline static double lastYValue = 0;
// 	inline static double lastZValue = 0;
// 	static void backup()
// 	{
// 		if (step == 2) logInfo("fields backup %g", lastXValue);
// 		lastXValue_copy = lastXValue;

// 	}
// 	static void restore()
// 	{
// 		if (step == 2) logInfo("fields restore was: %g rest: %g", lastXValue, lastXValue_copy);
// 		lastXValue = lastXValue_copy;
// 	}
// 	static double maxPulseDuration()
// 	{
// 		return XField::instance.maxPulseDuration();
// 	}
// 	//CALC
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		lastXValue = 0.0;
// 	}
// 	template <REP R, OPTIMS opt, typename F>
// 	double calc() const
// 	{
// 		return lastXValue;
// 	}
// 	template <REP R, typename WHEN> static constexpr inline bool canRun =
// 		COMPUTATION<double, FieldX>::template goodRep<R> && is_same_v<WHEN, LATE>;
// };

// template<typename ... Fields>
// struct LengthGauge : CouplingOperator<Fields...>
// {
// 	static constexpr REP couplesInRep = REP::X;
// 	static constexpr std::string_view CouplingName = "LG";
// 	static constexpr std::string_view name = "F_";
// 	using XField = typename CouplingOperator<Fields...>::XField;
// 	using YField = typename CouplingOperator<Fields...>::YField;
// 	using ZField = typename CouplingOperator<Fields...>::ZField;
// 	using CouplingOperator<Fields...>::lastXValue;
// 	using CouplingOperator<Fields...>::lastYValue;
// 	using CouplingOperator<Fields...>::lastZValue;

// 	// template <REP R, OPTIMS opt> inline double operator()(double x) const
// 	// {
// 	// 	if constexpr (is_same_v<XField, NoField>) return 0;
// 	// 	else if constexpr (R == REP::X) return PosOperator::template pos < R, opt >(i) * lastXValue;
// 	// 	else return 0.0;
// 	// }

// 	template <REP R, OPTIMS opt> inline double operator()(ind i, ind j = 0, ind k = 0) const
// 	{
// 		if constexpr (is_same_v<XField, NoField>) return 0;
// 		else if constexpr (R == REP::X)
// 		{
// 			if constexpr (is_same_v<YField, NoField>) return
// 				PosOperator::template pos < R, opt >(i) * lastXValue;
// 			else
// 			{
// 				if constexpr (is_same_v<ZField, NoField>) return
// 					(PosOperator::template pos < R, opt >(i) * lastXValue +
// 					 PosOperator::template pos < R, opt >(j) * lastYValue);
// 				else return (PosOperator::template pos < R, opt >(i) * lastXValue +
// 							 PosOperator::template pos < R, opt >(j) * lastYValue +
// 							 PosOperator::template pos < R, opt >(k) * lastZValue);
// 			}
// 		}
// 		else return 0.0;
// 	}
// 	template <REP R, OPTIMS opt>
// 	inline void precalc() const
// 	{
// 		if constexpr (R == REP::X)
// 		{
// 			if constexpr ((!is_same_v<XField, NoField>) && R == REP::X)
// 				lastXValue = XField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<YField, NoField>) && R == REP::X)
// 				lastXValue = YField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<ZField, NoField>) && R == REP::X)
// 				lastXValue = ZField::instance.operator()(timer);
// 		}
// 	}
// };

// template<typename ... Fields>
// struct CoulombGauge : CouplingOperator<Fields...>
// {
// private:
// 	inline static double lastTime_copy = 0;
// public:
// 	static constexpr REP couplesInRep = REP::P;
// 	static constexpr std::string_view CouplingName = "CG";
// 	static constexpr std::string_view name = "A_";
// 	using XField = typename CouplingOperator<Fields...>::XField;
// 	using YField = typename CouplingOperator<Fields...>::YField;
// 	using ZField = typename CouplingOperator<Fields...>::ZField;
// 	using CouplingOperator<Fields...>::lastXValue;
// 	using CouplingOperator<Fields...>::lastYValue;
// 	using CouplingOperator<Fields...>::lastZValue;
// 	inline static double lastTime = 0;
// 	static void backup()
// 	{
// 		CouplingOperator<Fields...>::backup();
// 		lastTime_copy = lastTime;
// 	}
// 	static void restore()
// 	{
// 		CouplingOperator<Fields...>::restore();
// 		lastTime = lastTime_copy;
// 	}
// 	// Full coupling: p^2/2 - p*A + A^2/2
// 	template <REP R, OPTIMS opt> inline double operator()(ind i, ind j = 0, ind k = 0) const
// 	{
// 		if constexpr (is_same_v<XField, NoField>) return 0;
// 		else if constexpr (R == REP::P)
// 		{
// 			if constexpr (is_same_v<YField, NoField>)
// 				return (-PosOperator::template pos < R, opt >(i) + 0.5 * lastXValue) * lastXValue;
// 			else
// 			{
// 				if constexpr (is_same_v<ZField, NoField>)
// 					return -((PosOperator::template pos < R, opt >(i) + 0.5 * lastXValue) * lastXValue +
// 							 (PosOperator::template pos < R, opt >(j) + 0.5 * lastYValue) * lastYValue);
// 				else return -((PosOperator::template pos < R, opt >(i) + 0.5 * lastXValue) * lastXValue +
// 							  (PosOperator::template pos < R, opt >(j) + 0.5 * lastYValue) * lastYValue +
// 							  (PosOperator::template pos < R, opt >(k) + 0.5 * lastZValue) * lastZValue);
// 			}
// 		}
// 		else return 0.0;
// 	}
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		//TODO: Check whether an even number of cycles was chosen to guarantee A=0 at the end
// 		CouplingOperator<Fields...>::template prepare<M, R, opt>();
// 		lastTime = 0.0;
// 	}
// 	// A(t) = -int_0^t E(t) dt
// 	template <REP R, OPTIMS opt>
// 	inline void precalc() const
// 	{
// 		if constexpr (R == REP::P)
// 		{
// 			if constexpr ((!is_same_v<XField, NoField>) && R == REP::P)
// 				lastXValue -= (timer - lastTime) * XField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<YField, NoField>) && R == REP::P)
// 				lastYValue -= (timer - lastTime) * YField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<ZField, NoField>) && R == REP::P)
// 				lastZValue -= (timer - lastTime) * ZField::instance.operator()(timer);
// 			lastTime = timer;
// 		}
// 	}
// };

// template<typename ... Fields>
// struct VelocityGauge : CouplingOperator<Fields...>
// {
// private:
// 	inline static double lastTime_copy = 0;
// public:
// 	static constexpr REP couplesInRep = REP::P;
// 	static constexpr std::string_view CouplingName = "VG";
// 	static constexpr std::string_view name = "A_";
// 	using XField = typename CouplingOperator<Fields...>::XField;
// 	using YField = typename CouplingOperator<Fields...>::YField;
// 	using ZField = typename CouplingOperator<Fields...>::ZField;
// 	using CouplingOperator<Fields...>::lastXValue;
// 	using CouplingOperator<Fields...>::lastYValue;
// 	using CouplingOperator<Fields...>::lastZValue;
// 	inline static double lastTime = 0;
// 	static void backup()
// 	{
// 		CouplingOperator<Fields...>::backup();
// 		lastTime_copy = lastTime;
// 	}
// 	static void restore()
// 	{
// 		CouplingOperator<Fields...>::restore();
// 		lastTime = lastTime_copy;
// 	}
// 	// Full coupling: p^2/2 - p*A + A^2/2
// 	template <REP R, OPTIMS opt> inline double operator()(ind i, ind j = 0, ind k = 0) const
// 	{
// 		if constexpr (is_same_v<XField, NoField>) return 0;
// 		else if constexpr (R == REP::P)
// 		{
// 			if constexpr (is_same_v<YField, NoField>)
// 				return (-PosOperator::template pos < R, opt >(i)) * lastXValue;
// 			else
// 			{
// 				if constexpr (is_same_v<ZField, NoField>)
// 					return -((PosOperator::template pos < R, opt >(i)) * lastXValue +
// 							 (PosOperator::template pos < R, opt >(j)) * lastYValue);
// 				else return -((PosOperator::template pos < R, opt >(i)) * lastXValue +
// 							  (PosOperator::template pos < R, opt >(j)) * lastYValue +
// 							  (PosOperator::template pos < R, opt >(k)) * lastZValue);
// 			}
// 		}
// 		else return 0.0;
// 	}
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		//TODO: Check whether an even number of cycles was chosen to guarantee A=0 at the end
// 		CouplingOperator<Fields...>::template prepare<M, R, opt>();
// 		lastTime = 0.0;
// 	}
// 	// A(t) = -int_0^t E(t) dt
// 	template <REP R, OPTIMS opt>
// 	inline void precalc() const
// 	{
// 		if constexpr (R == REP::P)
// 		{
// 			if constexpr ((!is_same_v<XField, NoField>) && R == REP::P)
// 				lastXValue -= (timer - lastTime) * XField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<YField, NoField>) && R == REP::P)
// 				lastYValue -= (timer - lastTime) * YField::instance.operator()(timer);
// 			if constexpr ((!is_same_v<ZField, NoField>) && R == REP::P)
// 				lastZValue -= (timer - lastTime) * ZField::instance.operator()(timer);
// 			lastTime = timer;
// 		}
// 	}
// };

// template<typename ... Fields>
// struct VelocityGaugeES : VelocityGauge<Fields...>
// {
// private:
// 	inline static double lastTime_copy = 0;
// public:
// 	static constexpr REP couplesInRep = REP::P;
// 	static constexpr std::string_view CouplingName = "VG_ES";
// 	static constexpr std::string_view name = "A_";
// 	using XField = typename CouplingOperator<Fields...>::XField;
// 	using YField = typename CouplingOperator<Fields...>::YField;
// 	using ZField = typename CouplingOperator<Fields...>::ZField;
// 	using CouplingOperator<Fields...>::lastXValue;
// 	using CouplingOperator<Fields...>::lastYValue;
// 	using CouplingOperator<Fields...>::lastZValue;
// 	// DIM==3: sqrt(2/3), DIM==2: sqrt(3)/2
// 	static constexpr double geoCoeff = DIM == 3 ? 0.81649658 : (DIM == 2 ? 0.8660254 : 1.0);
// 	inline static double lastTime = 0;

// 	constexpr VelocityGaugeES() :VelocityGauge<Fields...>() {}

// 	static void backup()
// 	{
// 		CouplingOperator<Fields...>::backup();
// 		lastTime_copy = lastTime;
// 	}
// 	static void restore()
// 	{
// 		CouplingOperator<Fields...>::restore();
// 		lastTime = lastTime_copy;
// 	}
// 	template <REP R, OPTIMS opt> inline double operator()(ind i, ind j = 0, ind k = 0) const
// 	{
// 		if constexpr (is_same_v<XField, NoField>) return 0;
// 		else if constexpr (R == REP::P)
// 		{
// 			if constexpr (is_same_v<YField, NoField>)
// 				return -((PosOperator::template pos < R, opt >(i)) * lastXValue +
// 						 (PosOperator::template pos < R, opt >(j)) * lastYValue +
// 						 (PosOperator::template pos < R, opt >(k)) * lastZValue);
// 			else
// 			{
// 				if constexpr (is_same_v<ZField, NoField>)
// 					return -((PosOperator::template pos < R, opt >(i)) * lastXValue +
// 							 (PosOperator::template pos < R, opt >(j)) * lastYValue);
// 				else return -((PosOperator::template pos < R, opt >(i)) * lastXValue +
// 							  (PosOperator::template pos < R, opt >(j)) * lastYValue +
// 							  (PosOperator::template pos < R, opt >(k)) * lastZValue);
// 			}
// 		}
// 		else return 0.0;
// 	}
// 	template <MODE M, REP R, OPTIMS opt>
// 	void prepare() const
// 	{
// 		//TODO: Check whether an even number of cycles was chosen to guarantee A=0 at the end
// 		CouplingOperator<Fields...>::template prepare<M, R, opt>();
// 		lastTime = 0.0;
// 	}
// 	// A(t) = -int_0^t E(t) dt
// 	template <REP R, OPTIMS opt>
// 	inline void precalc() const
// 	{
// 		if constexpr (R == REP::P)
// 		{
// 			if constexpr ((!is_same_v<YField, NoField>) && R == REP::P)
// 				lastYValue -= (timer - lastTime) * YField::instance.operator()(timer);
// 			else if constexpr ((!is_same_v<ZField, NoField>) && R == REP::P)
// 				lastZValue -= (timer - lastTime) * ZField::instance.operator()(timer);
// 			else if constexpr ((!is_same_v<XField, NoField>) && R == REP::P)
// 			{
// 				lastXValue -= geoCoeff * (timer - lastTime) * XField::instance.operator()(timer);
// 				lastYValue = lastXValue;
// 				lastZValue = lastXValue;
// 			}

// 			lastTime = timer;
// 		}
// 	}
// };

// template<typename ...Fields>
// struct LengthGauge_EckhardtSacha : CouplingOperator<Fields...>
// {
// 	static constexpr std::string_view CouplingName = "_LG_ES";
// 	static constexpr std::string_view name = "F_";
// 	using XField = typename CouplingOperator<Fields...>::XField;
// 	using YField = typename CouplingOperator<Fields...>::YField;
// 	using ZField = typename CouplingOperator<Fields...>::ZField;
// 	using CouplingOperator<Fields...>::lastXValue;
// 	using CouplingOperator<Fields...>::lastYValue;
// 	using CouplingOperator<Fields...>::lastZValue;
// 	double geoCoeff;

// 	constexpr LengthGauge_EckhardtSacha() :LengthGauge<Fields...>(), geoCoeff((ELEC == 2) ? sqrt(3.0) / 2.0 : sqrt(2.0 / 3.0)) {}

// 	template <REP R, OPTIMS opt> inline double operator()(ind i, ind j = 0, ind k = 0) const
// 	{
// 		if constexpr (is_same_v<XField, NoField>) return 0;
// 		else if constexpr (R == REP::X)
// 		{
// 			if constexpr (is_same_v<YField, NoField>) return PosOperator::template pos < R, opt >(i) * lastXValue;
// 			else
// 			{
// 				if constexpr (is_same_v<ZField, NoField>) return
// 					(PosOperator::template pos < R, opt >(i) * lastXValue +
// 					 PosOperator::template pos < R, opt >(j) * lastYValue);
// 				else return
// 					(PosOperator::template pos < R, opt >(i) * lastXValue +
// 					 PosOperator::template pos < R, opt >(j) * lastYValue +
// 					 PosOperator::template pos < R, opt >(k) * lastZValue);
// 			}
// 		}
// 		else return 0.0;
// 	}
// 	//PROP
// 	template <REP R, OPTIMS opt>
// 	inline void precalc() const
// 	{
// 		if constexpr (R == REP::X)
// 			if constexpr ((!is_same_v<XField, NoField>) && R == REP::X)
// 				lastXValue = geoCoeff * XField::instance.operator()(timer);
// 	}
// };








// EXCEPTIONAL PREDEFINED PULSE WITH ENVELOPE (for historical reasons)
// Use Field<ChemPhysEnvelope, ChemPhysPulse> instead
// double chp_tmp1, chp_tmp2, chp_tmp3;
// inline double ChemPhysPulseWithEnvelope(double time)
// {
// 	if (time < pulse_time)
// 	{
// 		chp_tmp1 = env_mult * timer;
// 		chp_tmp2 = omega * timer + phi - pi * ncycles;
// 		chp_tmp3 = sin(chp_tmp1);
// 		return F0 * chp_tmp3 * (chp_tmp3 * cos(chp_tmp2) + cos(chp_tmp1) * sin(chp_tmp2) / ncycles);
// 	}
// 	else return 0.0;
// }


// void initDefaultFieldVars()
// {
// 	phi0 = phi;
// 	phi = phi * 2.0 * pi;                       // "real" phase
// 	onecycle_time = 2.0 * pi / omega;           // au duration of one field cycle
// 	pulse_time = ncycles * onecycle_time; 	// au duration of the pulse
// 	env_mult = pi / pulse_time;	//used for Sin2Envelope //omega / (2.0 * ncycles)
// 	sin_pow = log(2) / log(1 / sin(pi * (ncycles - FWHM_percent) / ncycles / 2.0));
// 	gauss_pow = 4.0 * log(2) / POW2(FWHM_percent * onecycle_time);
// }