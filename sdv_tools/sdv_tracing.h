/*	Functions:	*/
//	initialize_tracing()-> Call 
//	define_regions("a","b",...) -> set event(1000,1)="a", event(1000,2)="b", ...
//	trace_event(x,y)		-> Add an event=x with value=y
//	enable_tracing()		-> Enable tracing from this point in time
//	disable_tracing()		-> Disable tracing from this point in time
/*-------------------*/


#ifdef EXTRAE
  #include "extrae_user_events.h"
  #define trace_event(x,y) Extrae_eventandcounters(x,y);
#elif defined(VEHAVE)
  #include <vehave-control.h>
  #define trace_event(x,y) vehave_trace(x,y);
#else
  #define trace_event(x,y);
#endif

#ifdef EXTRAE
	#define initialize_tracing() Extrae_init();
	#define define_regions(x,...)\
		{\
	  char * description_values[] = {"Zero", __VA_ARGS__};\
	  unsigned int nvalues=sizeof(description_values)/sizeof(char*);\
	  extrae_type_t type = x;\
	  extrae_value_t values[nvalues];\
		for(int i=0; i<nvalues; ++i) values[i]=i;\
	  Extrae_define_event_type (&type, "Code region", &nvalues, values, description_values);\
		}
	#ifdef EXTRAE_ALWAYS_TRACE
		#define enable_tracing();
		#define disable_tracing();
	#else
		#define enable_tracing() Extrae_restart();
		#define disable_tracing() Extrae_shutdown();
	#endif
#elif defined(VEHAVE)
	#define initialize_tracing();
	#define define_regions(x,...);
	#define enable_tracing() vehave_enable_tracing();
	#define disable_tracing() vehave_disable_tracing();
#else
	#define initialize_tracing();
	#define define_regions(x,...);
	#define enable_tracing();
	#define disable_tracing();
#endif
