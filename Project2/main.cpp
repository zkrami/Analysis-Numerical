#include "Interpolation.h"
using namespace Project2;
#pragma comment(linker, "/SUBSYSTEM:windows /ENTRY:mainCRTStartup")
int main(){

	Interpolation^ form = gcnew Interpolation(); 
	form->ShowDialog();

	return 0;
}