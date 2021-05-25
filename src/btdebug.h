#ifndef BT_DEBUG
#define BT_DEBUG

extern char debug_buffer[256];

void btdebug_init();
//void btdebug_print(const char *fmt, ...);
void btdebug_print(const char* msg);
void btdebug_enter(const char*);
void btdebug_exit(const char*);

double btdebug_ms();


#endif
