double* opt_kin = nullptr;		// array of e^(-ip^2/2) (constant!)
double* opt_vstat = nullptr;



// struct _Operator { OPTIMS maskOpt; };
template <REP R = REP::X>
struct IdentityOperator : _Operator
{
	static constexpr REP rep = REP::X;
	static constexpr string_view name = "ID";
	template <typename ... Args> auto operator()(Args...args) const noexcept { return 1.0; }
};

struct PosOperator : _Operator
{
	// double xmin;
	// double dx;
	// double pmin;
	// double dp;

	// PosOperator(double xmin, double dx, double pmin, double dp) :
	// 	xmin(xmin), dx(dx), pmin(pmin), dp(dp) {};

	// static constexpr string_view name = "Pos";
	// template <REP R, OPTIMS opt>  inline double pos(ind i)
	// {
	// 	if constexpr (opt & PRECOMP_COORDS)
	// 	{
	// 		if constexpr (bool(R & REP::X)) return xcoords[i];
	// 		if constexpr (bool(R & REP::P)) return pcoords[i];
	// 	}
	// 	else
	// 	{
	// 		if constexpr (bool(R & REP::X)) return xmin + i * dx;
	// 		if constexpr (bool(R & REP::P)) //WARNING if counting derivative! pos(i+1) -pos(i-1) 
	// 		{
	// 			if constexpr (bool(opt & UNTRANSPOSED_P_DIM))
	// 				return pmin + i * dp;
	// 			else
	// 			{
	// 				if (i >= n2) return dp * (i - n);//+ VelocityGauge<Field<1, SinPulse>>::lastXValue;
	// 				else return dp * i;// + VelocityGauge<Field<1, SinPulse>>::lastXValue;
	// 			}
	// 		}

	// 	}
	// }
};

struct KineticEnergyOperator : PosOperator
{
	// static constexpr REP rep = REP::P;
	// static constexpr string_view name = "KIN_EN";
	// constexpr KineticEnergyOperator() :PosOperator() {}

	// template <REP R, OPTIMS opt, typename ... Args>
	// double operator()(Args... args) const
	// {
	// 	if constexpr (opt & PRECOMP_KIN) return opt_kin[i];// + V_dynamic<R, opt>(pos<R, opt>(args)...);
	// 	else
	// 	{
	// 		return 0.5 * (POW2(pos<R, opt>(args)) + ...);// + V_dynamic<R, opt>(pos<R, opt>(args)...);
	// 	}
	// }
};

struct PotentialEnergyOperator : PosOperator
{
	// static constexpr REP rep = REP::X;
	// static constexpr string_view name = "POT_EN";
	// constexpr PotentialEnergyOperator() :PosOperator() {}

	// template < REP R, OPTIMS opt, typename ... Args>
	// double operator()(Args ... args) const
	// {
	// 	if constexpr (opt & PRECOMP_VSTAT) return opt_vstat[k + n * (i * n + j)];// + V_dynamic<R, opt>(pos<R, opt>(args)...);
	// 	else
	// 	{
	// 		return V_static(pos<R, opt>(args)...);// + V_dynamic<R, opt>(pos<R, opt>(args)...);
	// 	}
	// }
};

namespace Schrodinger
{
	/* Define named operators that can be used in computations */
	struct KinEnergy { };
	struct PotEnergy { };
	struct Couplings { };

	template <class V_Op, class C_Op, class GType>
	struct Spin0 : V_Op, C_Op, WF < Spin0<V_Op, C_Op, GType>, GType, 1>
	{
		using wf = WF < Spin0<V_Op, C_Op, GType>, GType, 1>;
		using wf::psi;
		using grid = typename wf::grid;
		using V = V_Op;
		using C = C_Op;
		Spin0() {}
		static constexpr REP couplesInRep = C::couplesInRep;
		// V pot;
		// static constexpr REP rep = REP::BOTH;
		// static constexpr string_view name = "Schrodinger";

		template < REP R, OPTIMS opt, typename ... Args>
		double operator()(Args ... args) const
		{
			if constexpr (R == REP::P)
			{
				if constexpr (couplesInRep == R)
					return (grid::kin_scale * (Power(grid::pos<R>(args), 2) + ...) +
							C_Op::template operator() < R, NO_OPTIMIZATIONS > (args...));
				else return grid::kin_scale * (Power(grid::pos<R>(args), 2) + ...);
			}
			else if (R == REP::X)
			{
				if constexpr (couplesInRep == R)
					return (V::template operator() < R, opt > (grid::pos<R>(args)...) +
							C_Op::template operator() < R, NO_OPTIMIZATIONS > (grid::pos<R>(args)...));
				else return V::template operator() < R, opt > (grid::pos<R>(args)...);
			}
		}
		template <MODE M>
		auto expOp(double val)
		{
			if constexpr (M == IM)
				return exp(-val);
			else
				return cos(-val) + I * sin(-val);
		}

		template <MODE M, REP R>
		void evolve(double delta)
		{
			for (i = 0; i < n0_l; i++)
			{
				if constexpr (DIM == 1)
					psi[i] *= expOp(delta * operator() < R, NO_OPTIMIZATIONS > (i + n0_o));
				else
				{
					readInd1 = i * grid::n;
					//Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to switch x<->y for DIM>1
					for (j = 0; j < grid::n; j++)
					{
						readInd2 = readInd1 + j;
						if constexpr (DIM == 2)
						{
							if (R == REP::X || MPI::region)
								psi[readInd2] *= te(delta * operator() < R, NO_OPTIMIZATIONS > (i + n0_o, j));
							else
								psi[readInd2] *= te(delta * operator() < R, NO_OPTIMIZATIONS > (j, i + n0_o));
						}
						else
						{
							readInd2 = readInd2 * grid::n;
							for (k = 0; k < grid::n; k++)
							{
								readInd3 = readInd2 + k;
								if (R == REP::X || MPI::region)
									psi[readInd3] *= te(delta * operator() < R, NO_OPTIMIZATIONS > (i + n0_o, j, k));
								else
									psi[readInd3] *= te(delta * operator() < R, NO_OPTIMIZATIONS > (j, i + n0_o, k));
							}
						}
					}
				}
			}
		}
	};
}
struct _XOperator : PosOperator { static constexpr string_view name = "X"; };
struct _YOperator : PosOperator { static constexpr string_view name = "Y"; };
struct _ZOperator : PosOperator { static constexpr string_view name = "Z"; };

template <typename V_Op>
struct PotentialDerivativeX :V_Op, _XOperator //TODO: HERE!!!
{
	static constexpr string_view name = "X_DerV";
	constexpr PotentialDerivativeX() : V_Op(), _XOperator() {}

	template <REP R, OPTIMS opt, typename ...Args>
	auto operator()(ind x, Args...args) const noexcept
	{
		return inv_2dx * (V_Op::template operator() < R, opt > (x + 1, args...) - V_Op::template operator() < R, opt > (x - 1, args...));
	}
};
template <typename V_Op>
struct PotentialDerivativeY : V_Op, _YOperator
{
	static constexpr string_view name = "Y_DerV";
	constexpr PotentialDerivativeY() : V_Op(), _YOperator() {}

	template <REP R, OPTIMS opt, typename ...Args>
	auto operator()(ind x, ind y, Args... args) const noexcept
	{
		return inv_2dx * (V_Op::template operator() < R, opt > (x, y + 1, args...) - V_Op::template operator() < R, opt > (x, y - 1, args...));
	}
};
template <typename V_Op>
struct PotentialDerivativeZ :V_Op, _ZOperator
{
	static constexpr string_view name = "Z_DerV";
	constexpr PotentialDerivativeZ() : V_Op(), _ZOperator() {}
	template <REP R, OPTIMS opt, typename ...Args>
	auto operator()(ind x, ind y, ind z) const noexcept
	{
		return inv_2dx * (V_Op::template operator() < R, opt > (x, y, z + 1) - V_Op::template operator() < R, opt > (x, y, z - 1));
	}
};



template <typename T_Op, typename V_Op, typename C_Op>
struct KineticOperator : T_Op, C_Op
{
	static constexpr string_view name = "T";
	double deltaMultiplier;
	static constexpr REP rep = REP::P;

	constexpr KineticOperator(double deltaMultiplier, double CAPlength) : T_Op(), deltaMultiplier(deltaMultiplier) {}
	template <MODE M, OPTIMS opt>

	auto operator()(ind i, ind j = 0, ind k = 0) const
	{
		double res = (T_Op::template operator() < REP::P, opt > (i, j, k) +
					  C_Op::template operator() < REP::P, opt > (i, j, k)) * dt * deltaMultiplier;
														  // if constexpr (opt & PRECOMP_KIN) return opt_kin[i];

		if constexpr (M == IM) return exp(-res);
		else return cos(-res) + I * sin(-res);
	}
	template <REP R, OPTIMS opt>
	inline void precalc() const
	{
		C_Op::template precalc<R, opt>();
	}
};

template <typename T_Op, typename V_Op, typename C_Op>
struct PotentialOperator : V_Op, C_Op
{
	static constexpr REP rep = REP::X;
	static constexpr string_view name = "V";

	double deltaMultiplier;

	constexpr PotentialOperator(double deltaMultiplier, double CAPlength) :
		V_Op(),
		C_Op(),
		deltaMultiplier(deltaMultiplier) {}

	template <MODE M, OPTIMS opt, typename ... Args>
	inline auto operator()(Args ... args) const
	{
		double res = (V_Op::template operator() < rep, opt > (args...) + C_Op::template operator() < REP::X, opt > (args...)) * dt * deltaMultiplier;
		// if constexpr (opt & PRECOMP_VSTAT) return opt_vstat[k + n * (i * n + j)];
		if constexpr (M == IM) return exp(-res);
		else return (cos(-res) + I * sin(-res));
	}
	template <REP R, OPTIMS opt>
	inline void precalc() const
	{
		C_Op::template precalc<R, opt>();
	}
};

template <typename T_Op, typename V_Op, typename C_Op>
struct PotentialOperatorWithCAP : PotentialOperator<T_Op, V_Op, C_Op>
{
	static constexpr string_view name = "U";
	inline static double CAPlength = 10.0;
	double eta;
	//The abs(deltaMultiplier) protects against blowup in certain split operator schemes
	constexpr PotentialOperatorWithCAP(double deltaMultiplier, double CAPlength) :
		PotentialOperator<T_Op, V_Op, C_Op>(deltaMultiplier, CAPlength),
		eta(CAPlength > 0 ? POW4(3.0 / CAPlength) : 0) {}

	// double operator ()(Args ... args) const noexcept;
	// template <typename ...Args>
	inline double cap_bare(double x, double y = 0, double z = 0) const noexcept
	{
		// double rdist, rdist2, rdist3;
		// double L2_effective = L / 2 - CAPlength;
		// rdist = (fabs(x) - L2_effective);
		// rdist2 = (fabs(y) - L2_effective);
		// rdist3 = (fabs(z) - L2_effective);
		// if (rdist3 <= 0 && rdist2 <= 0 && rdist <= 0) return 0.0;
		// else
		// {
		// 	rdist = rdist > 0 ? POW4(rdist) : 0;
		// 	rdist2 = rdist2 > 0 ? POW4(rdist2) : 0;
		// 	rdist3 = rdist3 > 0 ? POW4(rdist3) : 0;
		// 	return -(rdist3 + rdist2 + rdist) * eta;
		// }
		return 0.0;
	}
	template <MODE M, OPTIMS opt, typename ... Args>
	inline auto operator()(Args...args) const noexcept
	{
		// auto val = PotentialOperator<T_Op, V_Op, C_Op>::template operator() < M, opt, Args... > (args...);
		// if constexpr (bool(M & IM)) return val;
		// else return val * exp(cap_bare(PosOperator::pos<REP::X, opt>(args)...) * dt * PotentialOperator < T_Op, V_Op, C_Op>::deltaMultiplier);
	}
};





