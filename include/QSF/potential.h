struct _InteractionConfig
{
	/* If you think that a generic field parameter is lacking add it here */
	double Nsoft, Esoft, Ncharge, Echarge = -1;//,yukawaDecay;
};
struct InteractionPrototype : _InteractionConfig
{
	InteractionPrototype() = default;
	InteractionPrototype(_InteractionConfig p) : _InteractionConfig(p) {};
};

template <size_t ID, typename ... Inter>
struct Potential
{
	using type = Potential;
	inline static Potential instance;
	int count = sizeof...(Inter);

	InteractionPrototype* inters[sizeof...(Inter)]{};
	static constexpr std::string_view name = join_v<Inter::name...>;

	template <FREE_COORD fc, typename ... Args>
	inline double operator()(Args...coords)
	{
		// size_t i = 0;
		// double val = 0;
		// ([&] {
		// 	val += static_cast<Inter*>(inters[i])->operator() < fc > (coords...) + ...);
		// 	i++;
		// }(), ...);
		// return val;
		return (static_cast<Inter*>(inters[Index_v<Inter, Inter...>])->template operator() < fc > (coords...) + ...);
	}
};


struct CoulombInteraction : InteractionPrototype
{
	double NEcharge;
	static constexpr std::string_view name = "C";
	CoulombInteraction() = default;
	CoulombInteraction(_InteractionConfig p) :
		InteractionPrototype(p), NEcharge(Ncharge* Echarge) {}

	template <FREE_COORD fc, typename ... Args>
	inline double operator()(Args...coords)
	{
		//TODO: finish
		return NEcharge / sqrt(((coords * coords) + ...) + Nsoft);
	}
};

struct EckhardtSachaInteraction : InteractionPrototype
{
	double NEcharge;
	double EEcharge;
	static constexpr std::string_view name = "ES";
	EckhardtSachaInteraction() = default;
	EckhardtSachaInteraction(_InteractionConfig p) :
		InteractionPrototype(p), NEcharge(Ncharge* Echarge),
		EEcharge(Echarge* Echarge) {}

	inline double ee_EckhardtSacha(double x, double y)
	{
		return EEcharge / sqrt((x - y) * (x - y) + x * y + Esoft);
	}
	inline double Ne_EckhardtSacha(double x)
	{
		return NEcharge / sqrt(x * x + Nsoft);
	}

	template <FREE_COORD fc> inline double operator()(double x)
	{
		return Ne_EckhardtSacha(x);
	}
	template <FREE_COORD fc> inline double operator()(double x, double y)
	{
		if (fc == FREE_COORD::Y)
			return Ne_EckhardtSacha(x);
		if (fc == FREE_COORD::X)
			return Ne_EckhardtSacha(y);
		else
			return ee_EckhardtSacha(x, y) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y);
	}
	template <FREE_COORD fc>
	inline double operator()(double x, double y, double z)
	{
		if (fc == FREE_COORD::Z)
			return ee_EckhardtSacha(x, y) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y);
		else if (fc == FREE_COORD::Y)
			return ee_EckhardtSacha(x, z) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(z);
		else if (fc == FREE_COORD::X)
			return ee_EckhardtSacha(y, z) + Ne_EckhardtSacha(y) + Ne_EckhardtSacha(z);
		else if (fc == FREE_COORD::YZ)
			return Ne_EckhardtSacha(x);
		else if (fc == FREE_COORD::XZ)
			return Ne_EckhardtSacha(y);
		else if (fc == FREE_COORD::XY)
			return Ne_EckhardtSacha(z);
		else return ee_EckhardtSacha(x, y) + ee_EckhardtSacha(y, z) + ee_EckhardtSacha(z, x) + Ne_EckhardtSacha(x) + Ne_EckhardtSacha(y) + Ne_EckhardtSacha(z);
	}


};
