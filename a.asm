; Redkina Yana Course VO_Assembler Examination Project

COMMENT #
THE TASK: 
Solution of the Diophantine equation using genetic algorithm:
A1*x1 + A2*x2 + A3*x3  = D
xi - unknown positive integers (byts); Ai, D - given positive integer constants (bytes)
The size of the initial population N is set by the user in the range 4 <= N <= 10. The initial population is formed randomly.
Stop Criteria:
- exceeding the number of iterations M specified by the user;
- achieving a zero value of the objective function (residual equation).
Type of selection: proportional selection scheme [4]
Type of crossing: single-point [4] - solutions exchange one bit in any of xi
Mutation: change a randomly selected bit in any of xi. The probability of mutation is set by the user
#

include console.inc

;--------------------------------------VARIABLES-----------------------------------------;

.data

a = 3               ; const a for Lehmer algorithm
m = 251             ; const m for Lehmer algorithm
A1 = 1              ; const A1 of the Diophantine equation
A2 = 2              ; const A3 of the Diophantine equation
A3 = 3              ; const A4 of the Diophantine equation
D = 418            ; const D of the Diophantine equation
MaxN = 10           ; max number of population
X db MaxN dup (3 dup (?)) ; array of Xi variable
Y dw MaxN dup (?)   ; DELTA array for value of the objective function
U dw MaxN dup (?)   ; fitness array for value of the objective function
RS dw MaxN dup (?)  ; running sum with Fi/Fcp values
Rd db 30 dup (?)    ; array of random numbers
Z dw 1000,?,?,?     ; array for storing the best solution within iteration
N dd ?              ; input number of population
M dd ?              ; input number of max iteration
P dd ?              ; input number of mutation probability
PM db ?             ; input program mode
Sol dd 0            ; variable-mark if the objective function = zero then Sol:= 1


;--------------------------------------PROCEDURES----------------------------------------;

.code

;----- Procedure for random number generation -----;

LemherRNG proc

    ; prologue
    push ebp
    mov ebp, esp
    push ebx
    push ecx
    push eax
    push edx
    push esi
    mov ecx,[ebp+12]        ; first parameter (size of array)
    lea ecx,[ecx + ecx * 2] ; size * 3 - for 3 X
    mov ebx,[ebp+8]         ; second parameter (link to array)
    dec ecx                 ; ecx := N - 1, because first element of array already exits (the seed)
    
@1: xor edx, edx
    movzx eax, byte ptr [ebx]          ; X[i]
    add ebx, 1              ; i := i + 1
    mov si, a
    mul si                  ; X[i] * a   
    mov esi, m
    div esi                 ; X[i] * a mod m
    mov byte ptr[ebx], dl   ; X[i+1] := X[i] * a mod m
    loop @1

    ; epilogue
    pop esi
    pop edx
    pop eax 
    pop ecx
    pop ebx
    pop ebp
    ret 2*4

LemherRNG endp


;----- Procedure for calculating the value of the objective function -----;

ObjFunc proc

    ; prologue
    push ebp
    mov ebp, esp
    push eax
    push ebx
    push ecx
    push edx
    push esi
    push edi
    mov ecx,[ebp+12]        ; first parameter (X array with variables Xi)
    mov edx,[ebp+8]         ; second parameter (Y array with value of objective function)

    mov esi, 0              ; counter
@1:
    mov al, [ecx]
    mov bl, A1
    mul bl              ; X1*A1
    mov di, ax          ; result = X1*A1

    mov al, [ecx+1]
    mov bl, A2
    mul bl              ; X2*A2
    add di, ax          ; result = X1*A1 + X2*A2

    mov al, [ecx+2]
    mov bl, A3
    mul bl              ; X3*A3
    add di, ax          ; result = X1*A1 + X2*A2 + X3*A3

    mov ax, D
    sub di, ax          ; result = D - (X1*A1 + X2*A2 + X3*A3

    cmp di, 0           ; cheking if the result=0 ( then sol := 1)
    jne @2
    mov Sol, 1

@2: cmp di, 0
    jns @4
    neg di
@4:  mov word ptr [edx], di          ; storing value of function in Y array
;    outwordln word ptr [edx],, ' Delta= '  ; output intermediate calculations (value of objective function)

    cmp di, word ptr Z[0]           ; comparing new value of function with previos one
    jae @3
    
    mov word ptr Z[0], di           ; if new value below - then storing solution Xi in Z array
    xor ax, ax
    mov al, byte ptr [ecx]
    mov word ptr Z[2],  ax
    mov al, byte ptr [ecx+1]
    mov word ptr Z[4],  ax
    mov al, byte ptr [ecx+2]
    mov word ptr Z[6],  ax

@3: add ecx, 3                      ; next iteration
    add edx, 2
    inc esi
    cmp esi, N
    jne @1

    ; epilogue
    pop edi
    pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
    pop ebp
    ret 2*4

ObjFunc endp


;----- Procedure for output population -----;

Output proc

    ; prologue
    push ebp
    mov ebp, esp
    push ecx
    push edx
    push ebx
    mov ebx,[ebp+8]       ; first parameter (X array with population)

    outstrln ' '
    outstrln '     Population: '
    outstrln '--------------------'
    outstrln 'N    x1    x2    x3'
    outstrln '--------------------'
    mov ecx, N          ; counter
    mov edx, 0          ; index
@1: 
    outword edx         ; output number N and Xi
    outword byte ptr [ebx], 5, ' '
    outword byte ptr [ebx+1], 5, ' '
    outwordln byte ptr [ebx+2], 5, ' '

    add ebx, 3           ; next iteration
    inc edx
    dec ecx
    jnz @1
    outstrln '--------------------'

    ; epilogue
    pop ebx
    pop edx
    pop ecx
    pop ebp
    ret 1*4

Output endp

;----- Procedure "Selection Operator" -----;

SelectOp proc 

    ; prologue
    push ebp
    mov ebp, esp
    push eax
    push ebx
    push ecx
    push edx
    push esi
    push edi
    mov ecx,[ebp+12]        ; first parameter (Y array with values of the objective function Fi)

    ; calculating fitness (Fi/Fcp)

    xor ax, ax
    xor si,si
    mov edi, 0                  ; counter
@1: mov bx, word ptr [ecx] 
;    outwordln bx,, 'delta='
    xor dx, dx
    mov ax, 1000
    div bx                      ; forming value of function as inverse coefficients of fitness: 1000/Fi
    mov word ptr U[edi*2], ax   ; storing result in U array
    add si, ax                  ; forming total of all fitness coefficient : F=F1+F2...Fn
        
    add ecx, 2                  ; next iteration
    inc edi
    cmp edi, N
    jne @1

    xor ax, ax                  ; new loop
    mov ecx, 0                  ; index
    mov edi, N                  ; counter
@2: mov ax, word ptr U[ecx]
    mov bx, 100
    mul bx
    xor dx, dx
    mov bx, si
    div bx                      ; fitness coefficient for each solution n : Fi/Total
    cmp ax, 0
    jne @10
    inc ax
@10: mov word ptr U[ecx], ax     ; storing result in U array
;    outwordln word ptr U[ecx],, ' Fitness(%)='     ; output intermediate calculations (fitness coefficient)

    add ecx, 2                  ; next iteration
    dec edi
    jne @2

    ; calculating running sum

    mov ax, word ptr U[0]       ; first element in RS array
    mov word ptr RS[0], ax
;    outwordln word ptr RS[0],, ' Runnig Sum= '     ; output intermediate calculations (value of first running sum)

    xor ecx, ecx
    add ecx, 2                  ; index for U array with fitness indices
    mov edi, 2                  ; index for RS array with running sum
    mov esi, N                  ; counter
    dec esi                     ; N-1, because first element already exists
@3: mov ax, word ptr U[ecx]
    mov dx, word ptr RS[edi-2]
    add ax, dx                  ; result = RS[i]+U[i]
    mov word ptr RS[edi], ax    ; storing result in RS array
;    outwordln word ptr RS[edi],, ' Runnig Sum= '   ; output intermediate calculations (value of running sum)

    add edi, 2                  ; next iteration
    add ecx, 2
    dec esi
    jne @3

    ; Calling the Procedure LemherRNG :
    
@4: push 1                      ; first parameter to the stack
    push offset Rd              ; second parameter to the stack
    call LemherRNG
    
    mov esi,[ebp+8]             ; second parameter (Rd array with random numbers)
    mov al, byte ptr [esi+2]    ; change initial value for next call of random procedure
    cmp al, 0
    jne @11
    mov byte ptr [esi], 1
    jmp @12
@11:mov byte ptr [esi], al
     
;    outstr 'Random: '          ; output intermediate calculations (two random numbers)
;    outword byte ptr [esi+1],, ' '
;    outwordln byte ptr [esi+2],, ' '

@12: ; scaling random numbers
    
    mov ecx, 2  ; Counter 
    mov edi, N  ; calculating last index of RS
    lea edi, [edi+edi]
    mov bx, word ptr RS[edi-2]  ; runningSumMax
@13:movzx ax, byte ptr [esi+1]    ; random number for proc
    mul bx
    mov di, 251
    div di
    mov byte ptr [esi+1], al
    inc esi
    loop @13
     
;    mov esi,[ebp+8] 
;    outstr 'Random: '          ; output intermediate calculations (two random numbers)
;    outword byte ptr [esi+1],, ' '
;    outwordln byte ptr [esi+2],, ' '

    ; search random in range of running summ  

    mov edx,[ebp+8]         ; second parameter (Rd array with random numbers)
    inc edx                 ; taking the second element of the array with random numbers
    mov esi, N              ; counter for loop on RS array (0..N)
    mov edi, 2              ; counter for loop on Rd array (2)
    xor ecx, ecx            ; index for RS array (with running sum)
@9: movzx ebx, byte ptr [edx]   ; cheking out of range
    lea ecx, [esi-1]            ; index of last element in RS array
    movzx eax, word ptr RS[ecx*2]
    cmp ebx, eax
    jb @5
    mov byte ptr [edx], cl  ; storing index of last value of RS in Rd array

@5: xor ecx, ecx
@7: movzx eax, word ptr RS[ecx*2]
    cmp ebx, eax            ; finding right index in diaposone of RS array
    ja @6
    mov byte ptr [edx], cl  ; storing index as number for pair in Rd array
    jmp @8

@6: add ecx, 1              ; next iteration for RS array
    dec esi
    jne @7

@8: inc edx                 ; next iteration for Rd array
    xor ecx, ecx
    mov esi, N
    dec edi
    jne @9

    ;cheking if number i <> number j - then forming a pair

    mov edx,[ebp+8]
    mov al, byte ptr[edx+1]
    mov bl, byte ptr[edx+2]
    cmp al, bl              ; compare pair : n[i] and n[j]
    je @4                   ; if n[i] = n[j] then jump back to run the procedure again (for forming a new pair)

;    outstr 'Pair:  '       ; output intermediate calculations (selected pair for crosssing)
;    outword byte ptr [edx+1],, ' '
;    outwordln byte ptr [edx+2],, ' '
        
    ; epilogue 
    pop edi
    pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
    pop ebp
    ret 2*4

SelectOp endp


;----- Procedure "Crossing Operator" (in X3) -----;

CrossOp proc

    ; prologue
    push ebp
    mov ebp, esp
    push eax
    push ebx
    push ecx
    push edx
    push esi
    push edi  
    mov edi,[ebp+12]            ; first parameter (X array with variables Xi)
    mov esi,[ebp+8]             ; second parameter (Rd array with number of pair)

    mov al, byte ptr [esi+1]    ; calculating index of first element X3 in pair
    mov bl, 3
    mul bl
    add ax, 2                   
    mov byte ptr [esi+1], al    ; storing this index in Rd array
    lea edi, [edi+eax]
    movzx dx, byte ptr [edi]    ; value X3(i) with this index

    mov al, byte ptr [esi+2]    ; calculating index of second element X3 in pair
    mov bl, 3
    mul bl
    add ax, 2
    mov byte ptr [esi+2], al    ; storing this index in Rd array
    mov edi,[ebp+12]
    lea edi, [edi+eax]
    movzx cx, byte ptr [edi]    ; value X3(j) with this index

    ; calculating random bit n to exange between bytes of X3

    movzx ax, byte ptr [esi]    ; getting random number from array with random numbers
    mov bl, 17
    div bl
;    outwordln al,, 'random bit='   ; output intermediate calculations (random bit for crossing)
    movzx bx, al                ; n bit to exange

    ; crossing bits between bytes of X3(i) and X3(j)

    mov ax, dx                  ; X3(i)
    bt cx, bx                   ; bit definition in position n of X3(j)
    jc @1                       ; if bit = 1 then goto @1
    btr ax, bx                  ; if bit = 0 then set "0" in position n of X3(i)
    jmp @2
@1: bts ax, bx                  ; set "1" in position n of X3(i)

@2: movzx ebx, byte ptr [esi+1] ; calculating index of first element X3 in pair
    mov edi,[ebp+12]
    lea edi, [edi+ebx]
    movzx dx, byte ptr [edi]    ; copy of X3(i)
    mov byte ptr [edi], al      ; writing new value of X3(i)

    mov bl, byte ptr [esi+2]    ; calculating index of second element X3 in pair
    mov edi,[ebp+12]
    lea edi, [edi+ebx]
    movzx cx, byte ptr [edi]    ; X3(j)

    movzx ax, byte ptr [esi]    ; calculating same bit n to exange
    mov bl, 17
    div bl                      
    movzx bx, al                ; n bit to exange
    bt dx, bx                   ; bit definition in position n of X3(i)
    jc @3                       ; if bit = 1 then goto @3
    btr cx, bx                  ; if bit = 0 then set "0" in position n of X3(j)
    jmp @4
@3: bts cx, bx                  ; set "1" in position n of X3(j)
@4: mov byte ptr [edi], cl      ; writing new value of X3(j)
    
    ; epilogue
    pop edi
    pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
    pop ebp
    ret 2*4

CrossOp endp


;----- Procedure "Mutation Operator" (in X2) -----;

MutatOp proc

    ; prologue
    push ebp
    mov ebp, esp
    push eax
    push ebx
    push ecx
    push edx
    push esi
    push edi

    ; Calling the Procedure LemherRNG to get random number

    push N                  ; first parameter to the stack
    push offset Rd          ; second parameter to the stack
    call LemherRNG

    ; scaling random numbers

    mov ecx, N  ; Counter 
    lea ecx,[ecx + ecx * 2]
    dec ecx
    mov esi, 1
    mov bx, 100  ; random max
@5: movzx ax, byte ptr Rd[esi]    ; random number from proc
    mul bx
    mov di, 251
    div di
    mov byte ptr Rd[esi], al
    inc esi
    loop @5   

    ; comparing random number with entered probability
    
    mov edi,[ebp+8]                 ; second parameter (X array with variables Xi)
    inc edi                         ; shift pointer to X2 
    mov esi,[ebp+12]
    
    mov ecx, N                        ; counter for loop mutation for each solution N
    mov ebx, 1                        ; index in Rd array
@1: movzx eax, byte ptr Rd[ebx]
;    outwordln eax,, 'random number=' ; output intermediate calculations (random number for mutation)
    cmp eax, esi                      ; comparing values : if above - goto next solution
    ja @2

    ; calculating random bit for mutation

    movzx ax, byte ptr Rd[ebx+1]      ; getting random number from array with random numbers
    mov dl, 16
    div dl                            ; ah = number bit for mutation
;    outwordln al,, 'random bit='     ; output intermediate calculations (random bit for mutation)
    xor ah, ah

    ; mutation in random bit in X2

    movzx dx, byte ptr [edi]; X2
    bt dx, ax               ; bit definition in position n of X2
    jc @3                   ; if bit = 1 then goto @3
    bts dx, ax              ; if bit = 0 then set "1" in position n of X2
    jmp @4
@3: btr dx, ax              ; set "0" in position n of X2
@4: mov byte ptr [edi], dl  ; writing new value of X2
    
@2: add ebx, 2              ; next random number
    add edi, 3              ; next solution
    dec ecx
    jne @1

    ; epilogue
    pop edi
    pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
    pop ebp
    ret 2*4

MutatOp endp

;--------------------------------------START CALCULATION---------------------------------;

start:

;----- Input part of the code with cheking out of range -----;

    clrscr
    outstrln '==================================================='
    outstrln ' Diophantine equation to solve : x1+2*x2+3*x3 = 418' 
    outstrln '==================================================='

    outstr ' Input Size of Population N (from 4 to 10): '
    inint N
    mov eax, N
    cmp eax, 4
    js ErrInput
    cmp eax, 10
    ja ErrInput

    outstr ' Input Number of Iteration M (from 1 to 100): '
    inint M
    mov eax, M
    cmp eax, 1
    js ErrInput
    cmp eax, 100
    ja ErrInput

    mov edx, 0          ; counter for iterations

    outstr ' Input Probability of Mutation P (from 1 to 100 %): '
    inint P
    mov eax, P
    cmp eax, 1
    js ErrInput
    cmp eax, 100
    ja ErrInput

    outstr ' Input Program Mode (0-test mode/1-main mode): '
    inint PM
    mov al, PM
    cmp al, 0
    js ErrInput
    cmp al, 1
    ja ErrInput

;----- Calling Procedure "LemherRNG" for initial population generation -----;  

    mov byte ptr X[0], 19   ; initial value in X array (the seed for random algorithm)
    mov byte ptr Rd[0], 7  ; initial value in Rd array (the seed for random algorithm)
    push N                  ; first parameter to the stack
    push offset X           ; second parameter to the stack
    call LemherRNG


;----- Calling the Procedure "Output" if Program mode is TEST -----;

PrintRes:     
    mov al, PM              ; cheking Program mode
    cmp al, 0
    jne MainProgMod

    outstrln ' '
    outwordln edx,, ' Iteration N '
    push offset X           ; first parameter to the stack
    call Output


;----- Calling Procedure "ObjFunc" for calculating value of the objective function -----;

MainProgMod:
    push offset X       ; first parameter to the stack
    push offset Y       ; second parameter to the stack
    call ObjFunc


;----- Cheking stop criteria and displaing results -----;

    mov eax, Sol        ; cheking mark if Objective function = 0
    cmp eax, 0
    je IterCrit
    outstrln ' '
    outstrln '===================='
    outstrln '    Solution is: '
FinalRes:
    outwordln word ptr Z[2],, ' x1= '
    outwordln word ptr Z[4],, ' x2= '
    outwordln word ptr Z[6],, ' x3= '
    outwordln word ptr Z[0],, ' Objective Function = '
    outwordln edx,, ' Number of iteration: '
    exit

IterCrit:   
    cmp edx, M          ; cheking if number of iteration = max
    jb NxtIret
    outstrln ' '
    outstrln '===================='
    outstrln '  Best Solution is: '
    jmp FinalRes

;----- start of the next iteration -----;
 
NxtIret:

;----- Calling Procedure "SelectOp" to select a pair of solutions for further crossing -----;
    
    mov ecx, N
SlctLoop:
    push offset Y       ; first parameter to the stack
    push offset Rd      ; second parameter to the stack
    call SelectOp

;----- Calling Procedure "CrossOp" for crossing bits between pair of solutions -----;
    
    push offset X       ; first parameter to the stack
    push offset Rd      ; second parameter to the stack
    call CrossOp
    loop SlctLoop       ; loop N times to form N pairs

;---- Calling Procedure "MutatOp" for mutating bit and jumping back to estimate a new population ----;

    push P              ; first parameter to the stack
    push offset X       ; second parameter to the stack
    call MutatOp
    
    inc edx
    jmp PrintRes

;----- Exit from Program in case incorrect input -----;

ErrInput:
    outstrln 'Incorrect input. Try again'
    exit

end start  