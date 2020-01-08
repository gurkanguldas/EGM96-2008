
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class Main {
public static void main(String[] args) {
	String data ="E:\\Desktop\\DengelemedeOzelKonular\\Jeoit_Modelleme\\Odev_1\\data1.txt",
		    EGM ="E:\\Desktop\\DengelemedeOzelKonular\\Jeoit_Modelleme\\egm96_to360.ascii";
		    	
	Gravitation Jeoid = new Gravitation();
	Jeoid.setData(data);
	Jeoid.setEGM(EGM);
	double JeoPot[]=Jeoid.Potential(5);/**MAX:100**/
	
	DecimalFormat df = new DecimalFormat("#.#");
	df.setMinimumIntegerDigits(7);
	df.setMaximumFractionDigits(5);
	
	ElipsoidalGravitation Elipsoid = new ElipsoidalGravitation();
	Elipsoid.setData(data);
	double ElipPot[]=Elipsoid.Potential();
	

	
	for(int j=0 ; j<20 ; j++) 
		System.out.println("W["+(j+1)+"]=  "+df.format(JeoPot[j])+"     U["+(j+1)+"]=  "+df.format(ElipPot[j])+"     T["+(j+1)+"]=  "+(JeoPot[j]-ElipPot[j]));
		
}
}
