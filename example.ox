#include <oxstd.oxh>
#include <oxdraw.h>
#import <modelbase>
#include "prod_tools.ox"

main()
{
    decl asNamesConst, asNames, vDataConst=range(1,2)|range(3,4), vData;
    asNames = {"vX1","vX2"};
    asNamesConst = asNames;
    println("%r",{"Obs. 1", "Obs. 2"},
     "%c",asNamesConst,"%cf",{"%9.1f", "%9.1f"},vDataConst);
    vData=translog(vDataConst,0,&asNames);
    println("%r",{"Obs. 1", "Obs. 2"},
     "%c",asNames,"%cf",{"%9.1f"},
        vData);
    asNames=asNamesConst;
    vData=translog(vDataConst,-1,&asNames);
    println("%r",{"Obs. 1", "Obs. 2"},
     "%c",asNames,"%cf",{"%9.1f", "%9.1f"},
        vData);
    asNamesConst = {"vX1","vX2","vD1"};
    asNames=asNamesConst;
    vDataConst~=<0;1>;println(columns(vDataConst));
    println("%r",{"Obs. 1", "Obs. 2"},
     "%c",asNamesConst,"%cf",{"%9.1f", "%9.1f"},vDataConst);
    
    vData=translog(vDataConst,0,&asNames,range(0,1));
    println("%r",{"Obs. 1", "Obs. 2"},
     "%c",asNames,"%cf",{"%9.1f"},
        vData);

    println("db_expand is a helper function to save time in enumerating\n variables for use in modelbase-derived classes\n when no selection of lags is required.");
    println("%v",asNamesConst, " then becomes ", "%v",db_expand(asNamesConst));

/*             const editor = vscode.window.activeTextEditor;
            let cursorPosition = editor.selection.start;
            let wordRange = editor.document.getWordRangeAtPosition(cursorPosition);
            let highlight = editor.document.getText(wordRange);
            var oxdocFolder = OxBin_1.GetOxDocFolder();
            var url2 = ("file://"+oxdocFolder + "//oxstd.html#"+highlight);
            console.log(url2)
            open(url2);
 */


 }