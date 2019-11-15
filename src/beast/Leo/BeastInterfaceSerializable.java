/**
 * 
 */
/**
 * @author lr411
 *
 */
package beast.Leo;

import beast.core.BEASTInterface;

import java.io.FileOutputStream;
import java.io.ObjectOutputStream;

public interface BeastInterfaceSerializable extends BEASTInterface {
		public static void WriteObjectToFile(Object serObj, String filepath) {
			 
	        try {
	 
	            FileOutputStream fileOut = new FileOutputStream(filepath);
	            ObjectOutputStream objectOut = new ObjectOutputStream(fileOut);
	            objectOut.writeObject(serObj);
	            objectOut.close();
	 
	        } catch (Exception ex) {
	            ex.printStackTrace();
	        }
	    }
	   
	}
