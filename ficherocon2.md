<style>
.green {
    color: green;
    font-weight:700;
    font-size: 30px;
}
.heading1 {
    color: red;
    font-weight:700;
    font-size: 35px;
}
.heading2 {
    color: blue;
    font-weight:700;
    font-size: 30px;
}
border-collapse: collapse
</style>

<div class="green">
    Markdown css styles
</div>

# Esto es una prueba {#identifier .heading1}

Empieza cuando quieras, pero no te duermas escribiendo.

## Es un ejemplo {#identifier .heading2}

Pon ahi lo que tu quieras

|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table :is(td, th) {
                border: 1px solid black;
                padding: 0.3em;
                border-collapse: collapse
             }                
>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:16c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:45</td>  <td>0:02:35</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI108</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:28</td>  <td>0:01:22</td>  <td>0:00:53</td>  </tr>  </tbody>  </table>      | 0:10:20 |

|           |       hybrid (`mpirun` with "-sf" and "-pk")       |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  <th>h12r12t</th>  <th>h72r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:59</td>  <td>0:03:59</td>  <td>0:03:28</td>  <td>0:02:56</td>  <td>0:03:21</td>  <td>0:02:42</td>  <td>0:02:34</td>  <td>0:02:16</td> <td>0:01:29</td>  <td>0:04:02</td>  <td>0:01:34</td>  </tr>  </tbody>  </table>      |

|           |               hybrid (bare `mpirun`)               |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  <th>h12r12t</th>  <th>h72r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:40</td>  <td>0:03:58</td>  <td>0:02:23</td>  <td>0:02:10</td>  <td>0:02:03</td>  <td>0:02:23</td>  <td>0:01:49</td>  <td>0:01:47</td> <td>0:01:18</td>  <td>0:01:53</td>  <td>0:01:11</td>  </tr>  </tbody>  </table>      |

