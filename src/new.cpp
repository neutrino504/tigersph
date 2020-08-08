/*
 * new.cpp
 *
 *  Created on: Jul 18, 2020
 *      Author: dmarce1
 */
//#ifdef __AVX512F__
//#define ALIGN 64
//#else
//#define ALIGN 32
//#endif
//
//#include <cstdlib>
//#include <cstdio>
//
//void* operator new(std::size_t sz) {
//	void *ptr;
//	if ( posix_memalign(&ptr,ALIGN,sz) != 0) {
//		printf("Memory allcation error\n");
//		abort();
//	}
//	return ptr;
//
//}
//
//void operator delete(void* ptr) {
//	free(ptr);
//}
