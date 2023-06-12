# 普通物理學VPython作業

## 簡介

       這是某頂大電機系一年級必修，普通物理學的VPython作業。當你上網搜尋"Vpython GitHub"的時候，你應該已經知道是哪個大學了吧，而且你沒時間把它做完，這也是這裡存在的目的。簡單來說，這些檔案就是要送你學期成績12分;)。祝大學順利。

註: 以上檔案都有拿滿分，請放心"參考"。

## 目錄

###### 上學期

- HW1: Free Fall, Projectiles, Air Drags and Graph

- HW2: SHM, Collision, Pendulum and Newton Cradle

- HW3: Kepler's Planet Motion Law and Moon's Precession

- HW4: Phonons

- HW5: Equipartition of Energy of Molecular Gas
  
  包含兩個檔案: hw5 + hw5_diatomic

- HW6: Adiabatic Compression and Free Expansion of Ideal Gas
  
  包含兩個檔案: hw6 + hw6_histogram
  
###### 下學期

- HW7: Non-ideal Capacitor

- HW8: Drift Velocity

- HW9: Mutual Inductance

- HW10: RLC Circuit and Transient Response

- HW11: Imaging by a Thick Lens

- HW12: Diffraction from a Circular Aperture

## 使用建議

為了避免過多人同時"參考"，請記得稍加更改程式內容。例如:

1. 將變數名稱代換掉。特別注意程式的頭幾行通常是給定的數據，盡量不要更動這部分的變數名稱。

2. 增加註解

3. 將某些if的順序對調，例如:
   
   > if(a>0): A
   > 
   > else: B
   
   可以被改寫成
   
   > if(a<=0): B
   > 
   > else: A

4. 將**互相獨立**的程式順序對調，例如:
   
   > a = a + 2
   > 
   > b = b + 3
   
   可以被改寫成
   
   > b = b + 3
   > 
   > a = a + 2

5. 更改圖形的顏色、寬度

6. 增加無意義的空行

若不採取以上手段也無妨，因為待批改的作業數通常會多到無法檢查。

## Credits
這裡面大概一半都是抄前人的:) 我的工作是補足闕漏跟修正一些小bug、增加一些註解。
部分原作: https://github.com/rayray2002/no_plagiarism_vpython_sim

