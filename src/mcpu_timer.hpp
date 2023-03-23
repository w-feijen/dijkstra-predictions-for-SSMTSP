#include <boost/timer/timer.hpp>


class mcpu_timer { 	

	private: 

		// use cpu timer 
		boost::timer::cpu_timer timer;
	  	boost::timer::nanosecond_type total;
  		boost::timer::nanosecond_type last;

  		// count number of times a 0 time difference was measured
  		long int cnt_0;

	public:

  		mcpu_timer() : timer(), total(0), last(0), cnt_0(0) {}
  
  		// store current time
	  	void tik() { last = timer.elapsed().user; }  

	  	// add elapsed CPU time since last tik to total
  		void tok() { 
  			boost::timer::nanosecond_type diff;
  			diff = timer.elapsed().user - last;
  			if (diff == 0) cnt_0++;
  			total += diff;
  		}

  		// return number of times a 0 was measured 
  		long int zero_count() { return cnt_0; }

  		// returns total time measured in sec
		double out_sec() { return double(total) / (1000000000LL); }
		
  		// returns total time measured in nsec		
		boost::timer::nanosecond_type out_nsec() { return total; }

};


