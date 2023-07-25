|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:20c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:06</td>  <td>0:01:15</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:72c</th>  <th>pMPI:144c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:02:58</td>  <td>0:01:22</td>  <td>0:01:02</td>  </tr>  </tbody>  </table>      | 0:08:19 |

|          |                            hybrid                              |
|----------|----------------------------------------------------------------|
|t/min:sec | <table>  <thead>  <tr>  <th>h:18r/2t</th>  <th>h:9r/4t</th>  <th>h:6r/6t</th>  <th>h:12r/6t</th>  <th>h:18r/4t</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:04:31</th>  <th>0:04:12</th>  <th>0:06:05</th>  <th>0:05:22</th>  <th>0:04:13</th>  <th>0:02:48</td>  </tr>  </tbody>  </table>      |

|          |               Kokkos                     |
|----------|------------------------------------------|
|t/min:sec |  <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI:72c</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:06 </th>  <th>0:03:37</th>  <th>0:01:38</th>  <th>0:03:48</td>  </tr>  </tbody>  </table>      |                                     |
