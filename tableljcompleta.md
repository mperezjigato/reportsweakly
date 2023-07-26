|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:16c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:45</td>  <td>0:02:35</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI108</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:28</td>  <td>0:01:22</td>  <td>0:00:53</td>  </tr>  </tbody>  </table>      | 0:10:20 |

|           |       hybrid (`mpirun` with "-sf" and "-pk")       |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:59</td>  <td>0:03:59</td>  <td>0:03:28</td>  <td>0:02:56</td>  <td>0:03:21</td>  <td>0:02:42</td>  <td>0:02:34</td>  <td>0:02:16</td> <td>0:01:29</td>  </tr>  </tbody>  </table>      |

|           |               hybrid (bare `mpirun`)               |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:40</td>  <td>0:03:58</td>  <td>0:02:23</td>  <td>0:02:10</td>  <td>0:02:03</td>  <td>0:02:23</td>  <td>0:01:49</td>  <td>0:01:47</td> <td>0:01:18</td>  </tr>  </tbody>  </table>      |