/*
 this defines a code-block to be inserted in the neighbor search routines after the conditions for neighbor-validity are applied
 (valid particle types checked)
 */
  int numngb, no, p, task;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;
  MyDouble dx, dy, dz, dist;

#ifdef BOX_PERIODIC
  MyDouble xtmp;
#endif
#ifdef REDUCE_TREEWALK_BRANCHING
  t_vector box, hbox, vcenter;
#ifdef BOX_PERIODIC
  INIT_VECTOR3(boxSize_X, boxSize_Y, boxSize_Z, &box);
  INIT_VECTOR3(searchcenter[0], searchcenter[1], searchcenter[2], &vcenter);
  SCALE_VECTOR3(0.5, &box, &hbox);
#endif
#endif

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];
