#include <chrono>
#include <map>
#include <queue>
#include <stack>

namespace Timings
{
	using clock = std::chrono::high_resolution_clock;
	using TimeVar = clock::time_point;
	using std::chrono::seconds;
	using std::chrono::milliseconds;
	using std::chrono::duration;
	using std::chrono::duration_cast;
	TimeVar pass_start;

	constexpr int indentTotal = 36;
	constexpr int indentMult = 12;
	int indent = 0;
	class Timing
	{
	private:
		TimeVar start;
		duration<float, std::milli> elapsed = duration<float, std::milli>::zero();
		int indent = 0;
		bool globalTimer;
		bool passTimer;
		size_t ri;
		size_t pass;
	public:
		std::string_view name;
		// Main timer
		Timing() : name("MAIN"), globalTimer(true), passTimer(false)
		{
			start = clock::now();
		}
		// Pass timer
		Timing(std::string_view name, size_t ri, size_t pass, int indent) : name(name), globalTimer(false), passTimer(true), ri(ri), pass(pass), indent(indent* indentMult)
		{
			start = clock::now();

		}
		// Detail timer
		Timing(std::string_view name, int indent) : name(name), globalTimer(false), passTimer(false), indent(indent* indentMult)
		{
			// logInfo("CONS %s", name.data());
			start = clock::now();
		}
		void finish()
		{
			elapsed = clock::now() - start;
		}
		float rawTime() const
		{
			return elapsed.count();
		}
		void updateStart()
		{
			start = clock::now();
		}
		void additiveFinish()
		{
			// logInfo("CURRENT %s ms: %g", name.data(), elapsed.count());
			elapsed += clock::now() - start;
			// logWarning("CHANGED %s TO ms: %g", name.data(), elapsed.count());
		}
		void log()
		{
			if (globalTimer)
			{
				int ms = rawTime();
				int sec = ms / 1000;
				logSuccess("[%s] Total execution time: %4dh %3dmin %3dsec %3dms\n",
						   name.data(),
						   sec / 3600,
						   sec % 3600 / 60,
						   sec % 60,
						   ms % 1000);
			}
			else if (passTimer)
			{
				int ms = rawTime();
				int sec = ms / 1000;
				logSuccess("[ROUTINE %zu (%s) PASS %zu] Execution time: %4dh %3dmin %3dsec %3dms\n",
						   ri, name.data(), pass,
						   sec / 3600,
						   sec % 3600 / 60,
						   sec % 60,
						   ms % 1000);
			}
		}
		void logStats(float totalTime) const
		{
			if (passTimer)
			{
				int ms = rawTime();
				int sec = ms / 1000;
				logSuccess("ROUTINE %2zu PASS %2zu %-*s: %4dh %3dmin %3dsec %3dms (%10.8g %%)",
						   ri, pass, indentTotal - 19, name.data(),
						   sec / 3600,
						   sec % 3600 / 60,
						   sec % 60,
						   ms % 1000,
						   100.0 * rawTime() / totalTime);
			}
			else if (!globalTimer)
			{
				int ms = rawTime();
				int sec = ms / 1000;
				logSuccess("%*s%*s: %4dh %3dmin %3dsec %3dms (%10.5g %%)",
						   indent, name.data(), indentTotal - indent, "",
						   sec / 3600,
						   sec % 3600 / 60,
						   sec % 60,
						   ms % 1000,
						   100.0 * rawTime() / totalTime);
			}
		}
		// ~Timing() { logInfo("DEST %s", name.data()); }
	};

	std::stack<Timing> timings;
	std::queue<Timing> results;
	std::map<std::string_view, Timing> measures;
	void showTimings()
	{
		float totalTime = results.back().rawTime();
		logWarning("MEASUREMENTS: %s", results.back().name.data());
		for (auto iter = measures.rbegin(); iter != measures.rend(); ++iter)
			iter->second.logStats(totalTime);
		logWarning("TIME STATISTICS:");
		while (!results.empty())
		{
			results.front().logStats(totalTime);
			results.pop();
		}

	}
	void start(std::string_view name, size_t RI, size_t PASS)
	{
		pass_start = clock::now();
		timings.push(Timing(name, RI, PASS, indent));
		indent++;
	}
	namespace measure
	{
		void start(std::string_view name)
		{
			auto search = measures.find(name);
			if (search != measures.end()) search->second.updateStart();
			else measures.insert({ name, Timing(name, indent) });
		}
		void stop(std::string_view name)
		{
			auto search = measures.find(name);
			if (search != measures.end()) search->second.additiveFinish();
			else
			{
				logError("Measurement %s not found. Call Timings::measure::start(name) before ...stop(name).", name.data());
			}
		}
	}
	void start(std::string_view name)
	{
		timings.push(Timing(name, indent));
		indent++;
	}
	void start()
	{
		timings.push(Timing());
		// indent++;
	}

	void stop()
	{
		Timing& t = timings.top();
		t.finish();
		t.log();
		results.push(t);
		timings.pop();
		indent--;
		if (timings.empty()) showTimings();
	}
	int ETA(const size_t& step, const size_t& ntsteps)
	{
		if (step)
		{
			auto sec = duration_cast<seconds>(clock::now() - pass_start).count();
			if (step == 10)
				logInfo("ETA %g %td %zu %lld", (1.0 * ntsteps / step - 1.0) * sec, ntsteps, step, sec);
			return (1.0 * ntsteps / step - 1.0) * sec;
		}
		else return 0.0;
	}
}