function testVigraExtensions()
    numSuccess = 0;
    numTest = 0;
    
    disp('Testing: Loading/Creation of Input, Output');
    disp('         Input and Output Position (in Signature)');
    disp('         Typechecking of input');
    numSuccess = numSuccess + test_copy();
    numTest = numTest +1;
    
    disp('Testing: Loading/Creation of Input, Output');
    disp('         Input and Output Position (in Option Struct)');
    numSuccess = numSuccess + test_copy_option();
    numTest = numTest +1;
  
    disp('Testing: Loading/Creation of Input, Output');
    disp('         Input and Output Position (in Option Struct and Signature)');
    numSuccess = numSuccess + test_copy_mixed();
    numTest = numTest +1;
 
    
    
    
    disp('Testing: v_required');
    disp('         Argument Missing');
    numSuccess = numSuccess + test_v_required_in();
    numTest = numTest +1;
     
    disp('Testing: v_required');
    disp('         Argument Empty');
    numSuccess = numSuccess + test_v_required_in_empty();
    numTest = numTest +1;
    
    disp('Testing: v_required');
    disp('         Output Argument Missing');
    numSuccess = numSuccess + test_v_required_out();
    numTest = numTest +1;
    
    disp('Testing: v_optional');
    disp('         Input Argument');
    numSuccess = numSuccess + test_v_optional_in();
    numTest = numTest +1;
    
    disp('Testing: v_optional');
    disp('         Output Argument ');    
    numSuccess = numSuccess + test_v_optional_out();
    numTest = numTest +1;
    
    disp('Testing: v_default');
    disp('         Input Argument');
    numSuccess = numSuccess + test_v_default_in();  
    numTest = numTest +1;
    
    
    
    
    disp('Testing: ConstraintScalar');
    numSuccess = numSuccess + test_constr_scalar();
    numTest = numTest +1;
    
    disp('Testing: Enum');
    numSuccess = numSuccess + test_enum();
    numTest = numTest +1;
    
    disp('Testing: String');
    numSuccess = numSuccess + test_string();
    numTest = numTest +1;
    
    disp([num2str(numSuccess) ' of ' num2str(numTest) ' Tests successful.']);
return    

function success = test_copy_impl(handle)
    tMA     = handle(ones(3,3,3));
    tI      = handle(ones(3,3));
    tS      = handle(1);
    tV      = handle([1;1]);
    tC      = cell(4,1);
    success = true;
    
    [rMA, rI,rV rS, rS42, rC] = vigraTestCopy(tMA, tI,tV, tS, tC);
    if rMA(1) == handle(2)
        rMA(1) = 1;
    else
        rMA(1) = 2;
    end
    if rI(1) == handle(2)
        rI(1) = 1;
    else
        rI(1) = 2;
    end
    if rV(1) == handle(2)
        rV(1) = 1;
    else
        rV(1) = 2;
    end
    if rS(1) == handle(2)
        rS(1) = 1;
    else
        rS(1) = 2;
    end
    
    if ~strcmp(func2str(handle), class(rMA))
        success = false;
        disp(['...Output Class ' class(rMa) ' and  Input Class ' func2str(handle) ' do not match']);
    end
    
    if ~isequal(rMA, tMA);
        success = false;
        disp(['...Error in Loading/Saving MultiArray of Type ' class(tMA)]);
    end

    if ~isequal(rI, tI);
        success = false;
        disp(['...Error in Loading/Saving BasicImage of Type' class(tMA)]);
    end
    
    if ~isequal(rV, tV);
        success = false;
        disp(['...Error in Loading/Saving TinyVector of Type' class(tMA)]);
    end
        
    if rS ~= tS;
        success = false;
        disp(['...Error in Loading/Saving Scalar of Type' class(tMA) ' (Using pointers)']);
    end
    
    if rS42  ~= handle(42)
        success = false;
        disp(['...Error in Saving Scalar of Type' class(tMA)]);
    end
    
    if size(tC,1) ~= size(rC, 1) -1
        success = false;
        disp('...Error in Loading ConstCellArray');
    end
    

return

function success =  test_copy()
        success = false(10,1);
        success(1) = test_copy_impl(@double);
        success(2) = test_copy_impl(@single);
        success(3) = test_copy_impl(@uint8);
        success(4) = test_copy_impl(@uint16);
        success(5) = test_copy_impl(@uint32);
        success(6) = test_copy_impl(@uint64);
        success(7) = test_copy_impl(@int8);
        success(8) = test_copy_impl(@int16);
        success(9) = test_copy_impl(@int32);
        success(10) = test_copy_impl(@int64);
        success = isequal(success, true(10,1));
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end        
return


function success =  test_copy_option()
    tMA     = ones(3,3,3);
    tI      = ones(3,3);
    tS      = 1;
    tV      = [1;1];
    success = true;
    
    [rMA, rI,rV rS, rS42] = vigraTestCopyOpt(struct('field0',tMA, 'field1', tI,'field2', tV,'field3', tS));
    if rMA(1) == 2
        rMA(1) = 1;
    else
        rMA(1) = 2;
    end
    if rI(1) == 2
        rI(1) = 1;
    else
        rI(1) = 2;
    end
    if rV(1) == 2
        rV(1) = 1;
    else
        rV(1) = 2;
    end
    if rS(1) == 2
        rS(1) = 1;
    else
        rS(1) = 2;
    end

    
    if ~isequal(rMA, tMA);
        success = false;
        disp(['...Error in Loading/Saving MultiArray of Type ' class(tMA)]);
    end

    if ~isequal(rI, tI);
        success = false;
        disp(['...Error in Loading/Saving BasicImage of Type' class(tMA)]);
    end
    
    if ~isequal(rV, tV);
        success = false;
        disp(['...Error in Loading/Saving TinyVector of Type' class(tMA)]);
    end
        
    if rS ~= tS;
        success = false;
        disp(['...Error in Loading/Saving Scalar of Type' class(tMA) ' (Using pointers)']);
    end
    
    if rS42  ~= 42
        success = false;
        disp(['...Error in Saving Scalar of Type' class(tMA)]);
    end

    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');  
return

function success =  test_copy_mixed()
    tMA     = ones(3,3,3);
    tI      = ones(3,3);
    tS      = 1;
    tV      = [1;1];
    success = true;
    
    [rMA, rI,rV rS, rS42] = vigraTestCopyMixed(tMA, tI, struct('field2', tV,'field3', tS));
    if rMA(1) == 2
        rMA(1) = 1;
    else
        rMA(1) = 2;
    end
    if rI(1) == 2
        rI(1) = 1;
    else
        rI(1) = 2;
    end
    if rV(1) == 2
        rV(1) = 1;
    else
        rV(1) = 2;
    end
    if rS(1) == 2
        rS(1) = 1;
    else
        rS(1) = 2;
    end

    
    if ~isequal(rMA, tMA);
        success = false;
        disp(['...Error in Loading/Saving MultiArray of Type ' class(tMA)]);
    end

    if ~isequal(rI, tI);
        success = false;
        disp(['...Error in Loading/Saving BasicImage of Type' class(tMA)]);
    end
    
    if ~isequal(rV, tV);
        success = false;
        disp(['...Error in Loading/Saving TinyVector of Type' class(tMA)]);
    end
        
    if rS ~= tS;
        success = false;
        disp(['...Error in Loading/Saving Scalar of Type' class(tMA) ' (Using pointers)']);
    end
    
    if rS42  ~= 42
        success = false;
        disp(['...Error in Saving Scalar of Type' class(tMA)]);
    end
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');    
return


function success = test_v_required_in()
    tMA     = handle(ones(3,3,3));
    tI      = handle(ones(3,3));
    tS      = handle(1);
    tV      = handle([1;1]);
    success = false;
    try
        [rMA, rI,rV rS, rS42, rC] = vigraTestCopy(tMA, tI,tV, tS);
        disp('...Error - v_required not throwing exception if argument not given');
    catch
        success = true;
    end
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end   
    disp(' ');    
return;
   

function success = test_v_required_in_empty()
    tMA     = handle(ones(3,3,3));
    tI      = handle(ones(3,3));
    tS      = handle(1);
    tV      = handle([1;1]);
    tC      = cell(4,1);
    success = false;
    try    
        [rMA, rI,rV rS, rS42, rC] = vigraTestCopy(tMA, tI,tV, [], tC);
        disp('...Error - v_required not throwing exception if in argument empty');
    catch
        success = true;
    end
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end    
    disp(' ');    
return;


function success = test_v_required_out()
    tMA     = handle(ones(3,3,3));
    tI      = handle(ones(3,3));
    tS      = handle(1);
    tV      = handle([1;1]);
    tC      = cell(4,1);
    success = false;
    try       
        [rMA, rI,rV rS, rS42] = vigraTestCopy(tMA, tI,tV, tS);
        disp('...Error - v_required not throwing exception if out argument not given');
    catch
        success = true;
    end
    
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end    
    disp(' ');
return;
    
function success = test_enum()

    success = true;
    f = vigraTestEnum('first');
    if f ~= 1
        success = false;
        disp('...Error - Enumvalue wrong');       
    end
        
    s = vigraTestEnum('second');
    if s ~= 2
        success = false;
        disp('...Error - Enumvalue wrong');       
    end
    
    t = vigraTestEnum('third');
    if t ~= 3
        success = false;
        disp('...Error - Enumvalue wrong');       
    end
    if success == true
        success = false;
        try

            t = vigraTestEnum('forth');
            disp('...Error - no Error thrown if Enum does not exist');
        catch
            success = true;
        end
    end
 
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');
return;
    
function success = test_constr_scalar()
    success = false(11, 1);
    %ok
    try
        vigraTestConstrScalar(1, 3, 2, 2);
        success(1) = true;
    catch
        disp('...Error - Exception thrown with valid input')
    end
    
    %ok
    try
        vigraTestConstrScalar(false, 3,2,2);
        success(2) = true;
    catch
        disp('...Error  Exception thrown with valid logical input-')
    end
    
    %notok
    try
        vigraTestConstrScalar(3, 3,2,2);
        disp('...Error - No Exception thrown with invalid logical input')    
    catch
        success(3) = true;
    end    
    
    
    %notok
    try
        vigraTestConstrScalar(1, 4, 2, 2);
        disp('...Error - No Exception thrown with MinMax Constrained Scalar (Max exceeded)')
    catch
        success(4) = true;
    end    
    
    %notok
    try
        vigraTestConstrScalar(1, 1, 2, 2);    
        disp('...Error - No Exception thrown with MinMax Constrained Scalar (Min exceeded)')    
    catch
        success(5) = true;
    end    

    
    %ok
    try
        vigraTestConstrScalar(1, 3, 4, 2);
        success(6) = true;
    catch
        disp('...Error - Exception thrown with Valid Vals Constrained Scalar')
    end
    
    %notok
    try
        vigraTestConstrScalar(1, 3, 5, 2);
        disp('...Error - No Exception raised with invalid Vals Constrained Scalar')    
    catch
        success(7) = true;
    end    

    
    %ok
    try    
        vigraTestConstrScalar(1, 3, 2, 4);
        success(8) = true;
    catch
        disp('...Error - Exception thrown with Valid 2D3D Vals Constrained Scalar (3D)')
    end
    
    %ok
    try        
        vigraTestConstrScalar(1, 2, 2, 3);
        success(9) = true;
    catch
        disp('...Error - Exception thrown with Valid 2D3D Vals Constrained Scalar (2D)')
    end
    
    %ok
    try        
        vigraTestConstrScalar(1, 2, 2, 5);
        success(10) = true;
    catch
        disp('...Error - Exception thrown with Valid 2D3D Vals Constrained Scalar (2D)')
    end
    
    %notok
    try
        vigraTestConstrScalar(1, 2, 2, 4);
        disp('...Error - No Exception raised with invalid 2D3D Vals Constrained Scalar (2D)')    
    catch
        success(11) = true;
    end    

    %notok
    try
        vigraTestConstrScalar(true, 3, 2, 5);
        disp('...Error - No Exception raised with invalid 2D3D Vals Constrained Scalar (3D)')    
    catch
        success(12) = true;
    end       
    
    success = isequal(success, true(12, 1));

    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');    
return


function success = test_v_default_in()
    [a b c] = vigraTestDefault();
    success = true;
    if a ~= 2
        success = false;
        disp('...Error - Wrong default value for Scalar input');
    end
    
    if b ~= 2
        success = false;
        disp('...Error - Wrong 2D 3D default value (2D)');
    end
    
    if ~isequal(zeros(3,3), c)
        success = false;
        disp('...Error - Wrong default Array');
    end
    
    [a b c] = vigraTestDefault(3);
    if a ~= 3
        success = false;
        disp('...Error - Wrong set value for first input');
    end
    
    if b ~= 3
        success = false;
        disp('...Error - Wrong 2D 3D default value (3D)');
    end
    
    if ~isequal(zeros(3,3), c)
        success = false;
        disp('...Error - Wrong default Array');
    end
    
    [a b c] = vigraTestDefault(3, 7, ones(3,2));
    if a ~= 3
        success = false;
        disp('...Error - Wrong set value for first input');
    end
    
    if b ~= 7
        success = false;
        disp('...Error - Wrong 2D 3D default value (3D)');
    end
    
    if ~isequal(ones(3,2), c)
        success = false;
        disp('...Error - Wrong default Array');
    end    


    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');    
return


function success = test_v_optional_in()
    success = false(2,1);
    try
        [a b]  = vigraTestOptional(3, ones(3,3));
        success(1) = true;
        if a~= 3
            success(2) = false;
            disp('...Error - Error wrong value for Scalar input when Optional Parameter not supplied')
        end
        
        if ~isequal(b, ones(3,3))
            success(2) = false;
            disp('...Error - Error wrong value for Array input when Optional Parameter not supplied')
        end                
    catch
        disp('...Error - Exception Thrown when Optional Parameter supplied');
    end
    
    try
        [a b] = vigraTestOptional();
        success(2) = true;
        if a~= 2
            success(2) = false;
            disp('...Error - Error wrong value for Scalar input when Optional Parameter not supplied')
        end
        
        if ~isequal(b, zeros(3,3))
            success(2) = false;
            disp('...Error - Error wrong value for Array input when Optional Parameter not supplied')
        end        
    catch
        disp('...Error - Exception thrown when Optional Parameter not supplied');
    end
    
    success = isequal(success, true(2,1));
 

    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');     
return


function success = test_v_optional_out()
    success = false(2,1);
    try
        [a b]  = vigraTestOptional(3, ones(3,3));
        success(1) = true;
        if a~= 3
            success(2) = false;
            disp('...Error - Error wrong value for Scalar input when Optional Parameter not supplied')
        end
        
        if ~isequal(b, ones(3,3))
            success(2) = false;
            disp('...Error - Error wrong value for Array input when Optional Parameter not supplied')
        end                
    catch
        disp('...Error - Exception Thrown when Optional Output Parameter supplied');
    end
    
    try
        vigraTestOptional();
        success(2) = true;    
    catch
        disp('...Error - Exception thrown when Optional Output Parameter not supplied');
    end
    
    success = isequal(success, true(2,1));
 

    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end
    disp(' ');     
return


function success = test_string()
    success = true;
    a = vigraTestString();
    b = vigraTestString('user');
    
    if a ~= 1
        disp('...Error -default String not loaded');
        success = false;
    end

        
    if b ~= 2
        disp('...Error -user supplied string not loaded');
        success = false;
    end
 
    if success
        disp('                                                             Success');
    else
        disp('                                                             Fail');
    end    
    disp(' ');
return