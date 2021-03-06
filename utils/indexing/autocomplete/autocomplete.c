/*
 * Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* autocomplete.c -- generates autocomplete suggestion dictionary. In C
 * for speed as it traverses all indices.
 *
 * link with libm (math) ie -lm
 *
 * Author: Dan Shepaprd (ds23)
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>

#define READBUFFER 65536

/* Fields to exclude from autocomplete altogether (with 0 on end of list) */
char * bad_fields[] = {"domain_url",0};

/* Checks for membership of list of strings. If present reaturns 0, else
 * returns new length. Compact representation is \0 terminated strings,
 * followed by extra \0. NULL is acceptable as zero-length list. Cannot
 * store empty string.
 */

int strl_member(char **str,char *data,int max) {
  char *c,*d;
  int len,n=0;
  ptrdiff_t p;

  len = strlen(data);
  if(!len)
    return 0;
  if(!*str && max>0) { /* Initial malloc of zero-length list */
    *str = malloc(1);
    **str = '\0';
  }
  c=*str;
  while(*c) {
    for(d=data;*c && *c == *d;c++,d++)
      ;
    if(!*c && !*d)
      return 0;
    for(;*c;c++)
      ;
    c++;
    n++;
  }
  if(n<max) {
    p = c-*str; /* Remap c after realloc, also = length of alloc-1 */
    *str = realloc(*str,p+len+2);
    c = *str+p;
    strcpy(c,data);
    *(c+len+1) = '\0';
    n++;
  }
  return n;
}

int max_index_size = 100000;

int in_good_field = 0;
int good_field(char *name) {
  char **b;

  for(b=bad_fields;*b;b++)
    if(!strcmp(name,*b))
      return 0;
  return 1;
}

int quick_hash(char *str,int seed,int mod) {
  unsigned int hash = 0;

  for(;*str;str++)
    hash = (hash * 37 + (*str ^ seed));
  return hash % mod;
}

void process_tag(char *data) {
  char * field,*f;
  int idx;

  if(!strncmp(data,"field ",6)) {
    /* FIELD */
    if((field = strstr(data,"name=\""))) {
      f = index(field+6,'"');
      if(f)
        *f = '\0';
      if(good_field(field+6))
        in_good_field = 1;
    }
  } else if(!strcmp(data,"/field")) {
    in_good_field = 0;
  }
}

#define EFFORTLIMIT 12
#define MINTAIL 4

#define MAXMULT 4

char * make_raw_effort(char *data) {
  char *effort,*c,*d;
  int i,num;

  effort = malloc(strlen(data)*MAXMULT+1);
  for(c=effort,d=data;*d;d++) {
    num = 1;
    if(isdigit(*d))
      num = 4;
    for(i=0;i<num;i++)
      *(c++) = *d;
  }
  *c = '\0';
  return effort;
}

char * make_effort(char *data) {
  char *effort,*c,*d;
  int limit,i,num=0,len,off;

  limit = EFFORTLIMIT;
  len = strlen(data);
  /*
  if(limit > len)
    limit = len - MINTAIL;
  if(limit < 0)
    limit = 0;
  */
  off = len;
  for(i=0;i<len;i++) {
    if(isdigit(data[i]))
      num += 4;
    else
      num++;
    if(num>limit) {
      off = i+1;
      break;
    }
  }
  if(off<1) off=1;
  effort = malloc(off+1);
  strncpy(effort,data,off);
  effort[off]='\0';
  return effort;
}

int freq[] = {
  25,42,36,31,21,37,38,28,26,68,49,32,36,
  26,25,40,67,28,27,23,35,45,38,63,38,72
};

int annoyingness(char *data) {
  int n=0,len=0,v=0,f=0,letlen=0;
  char *c;

  for(c=data;*c;c++) {
    len++;
    if(!isalpha(*c))
      n+=100;
    if(strspn(c,"aeiou") > 2)
      n+=50; /* Too many vowels in a row */
    if(strcspn(c,"aeiou") > 3)
      n+=10; /* Too many consonants in a row */
    if(strspn(c,"aeiou"))
      v=1;
    if(*c>='a' && *c<='z') {
      f += freq[*c-'a'];
      letlen++;
    }
  }
  if(letlen) f/= letlen; else f = 100;
  if(f>30)
    n += (f-30)*5; /* unusual letters */
  if(!v)
    n += 30; /* no vowels */

  if(!len) return 0;
  return n/len;
}

#define MAX_ANNOYINGNESS 200

int annoyingness_size[MAX_ANNOYINGNESS];

void reset_annoyingness() {
  int i;

  for(i=0;i<MAX_ANNOYINGNESS;i++)
    annoyingness_size[i] = 0;
}

void register_annoyingness(int ann) {
  int i;

  for(i=ann;i<MAX_ANNOYINGNESS;i++)
    annoyingness_size[i]++;
}

int last_val = -1;
int annoyingness_threshold(int num) {
  int i;

  /* This method is on the critical path, so use a cache */
  if(last_val>=0) {
    if(last_val == MAX_ANNOYINGNESS-1 && annoyingness_size[MAX_ANNOYINGNESS-1]<num)
      return MAX_ANNOYINGNESS-1;
    if(annoyingness_size[last_val]<=num && annoyingness_size[last_val+1]>num)
      return last_val;
  }

  for(i=1;i<MAX_ANNOYINGNESS;i++) {
    if(annoyingness_size[i]>num) {
      last_val = i-1;
      return i-1;
    }
  }
  last_val = MAX_ANNOYINGNESS-1;
  return MAX_ANNOYINGNESS-1;
}

#define NUMEL 4
#define NUMFP 4 
struct counter {
  char * prefix;
  char *words;
  struct counter * next;
};

struct prefix_counter {
  int size,num;
  struct counter ** counter;
};

struct prefix_counter *prefixes,*words;

struct prefix_counter * make_prefix_counter(int size) {
  struct prefix_counter *out;
  int i;

  out = malloc(sizeof(struct prefix_counter));
  if(size) {
    out->size = size;
    out->num = 0;
    out->counter = malloc(sizeof(struct counter *)*size);
    for(i=0;i<size;i++)
      out->counter[i] = 0; 
  } else {
    *out = (struct prefix_counter){0,0,0};
  }
  return out;
}

void boost_prefix_counter(struct prefix_counter *in) {
  struct prefix_counter *out;
  struct counter *c,*d;
  int i,hash;

  out = make_prefix_counter(in->size*3/2+16);
  for(i=0;i<in->size;i++)
    for(c=in->counter[i];c;c=d) {
      d = c->next;
      hash = quick_hash(c->prefix,0,out->size);
      c->next = out->counter[hash];
      out->counter[hash] = c;
    }
  out->num = in->num;
  *in = *out;
}

long long int stat_naughty=0,stat_good=0,stat_words=0;
#define NAUGHTY_THRESHOLD 12
/* 1 = new, 0 = old/naughty */
int inc_counter(struct prefix_counter *pc,char *prefix,char *word) {
  int hash,whash,i,myloc,j,any,num;
  struct counter *c,*rec=0;

  if(pc->num >= pc->size/3)
    boost_prefix_counter(pc);
  hash = quick_hash(prefix,0,pc->size);
  for(c=pc->counter[hash];c;c=c->next) {
    if(!strcmp(c->prefix,prefix))
      rec = c;
  }
  if(!rec) {
    c = malloc(sizeof(struct counter));
    c->prefix = malloc(strlen(prefix)+1);
    strcpy(c->prefix,prefix);
    c->words = 0;
    c->next = pc->counter[hash];
    pc->counter[hash] = c;
    pc->num++;
    rec = c;
  } else if(!rec->words) {
    stat_naughty++;
    return 0;
  }
  num = strl_member(&(rec->words),word,NAUGHTY_THRESHOLD);
  if(num>=NAUGHTY_THRESHOLD) {
    free(rec->words);
    rec->words = 0;
    stat_naughty++;
    return 0;
  }
  stat_good++;
  if(num==0) {
    return 0;
  }
  return 1;
}

double stat_chlen = 0.0;

void dump_words(struct prefix_counter *pc) {
  struct counter *c;
  int i,thresh,ann;
  double value; 
  char *d;

  thresh = annoyingness_threshold(max_index_size);
  for(i=0;i<pc->size;i++)
    for(c=pc->counter[i];c;c=c->next) {
      d=c->words;
      if(d) {
        while(*d) {
          ann = annoyingness(d);
          if(ann<thresh) {
            value = ann*strlen(d);
            if(value > 0.5) {
              value = 12.0 - log10(value);
            } else {
              value = 12.0;
            }
            if(strlen(d)>6)
              value -= 0.1 * (strlen(d)-6);
            printf("%s\t%1.1f\n",d,value);
          }
          d += strlen(d)+1;
        }
      }
    }
}

void process_word(char *data) {
  char *effort;
  int thresh,ann;

  if(!*data)
    return;
  stat_words++;
  effort = make_effort(data);
  if(inc_counter(prefixes,effort,data)) {
    thresh = annoyingness_threshold(max_index_size);
    ann = annoyingness(data);
    if(ann<thresh) {
      register_annoyingness(ann);
    }
  }
  free(effort);
}

int isseparator(char c) {
  return c == '/' || c == ':' || c == ';' || c == '-' || c == '_' || c == '(' || c == ')';
}

void process_text(char *data) {
  int i;
  char *c,*d;

  if(!in_good_field)
    return;
  /* Remove punctuation attached to starts and ends of words */
  d = data;
  for(c=data;*c;c++) {
    if(ispunct(*c)) {
      if(c==data || !*(c+1) ||
         !isalnum(*(c-1)) || !isalnum(*(c+1)) ||
         isseparator(*c))
        *c = ' ';
    }
    if(isspace(*c)) {
      *c = '\0';
      if(*d)
        process_word(d);
      d = c+1;
    } else {
      *c = tolower(*c);
    }
  }
  process_word(d);
}

int tag_mode=0;
char *tagstr = 0;

void lex_part(char *part,int tag) {
  struct strings *s;

  if(tag_mode != tag) {
    /* Do stuff */
    if(tag_mode) {
      process_tag(tagstr);
    } else {
      process_text(tagstr);
    }
    free(tagstr);
    tagstr = 0;
    tag_mode = tag;
  }
  if(!tagstr) {
    tagstr = malloc(1);
    tagstr[0] = '\0';
  }
  if(strlen(part)) {
    tagstr = realloc(tagstr,strlen(tagstr)+strlen(part)+1);
    strcat(tagstr,part);
  }
}

int in_tag = 0;
/* at top level we just extract tag / non-tag and pass it down */
void lex(char *data) {
  int more,i;
  char match,*hit;

  while(*data) {
    hit = index(data,in_tag?'>':'<');
    if(hit)
      *hit = '\0';
    lex_part(data,in_tag);
    if(hit) {
      in_tag = !in_tag;
      data = hit+1;
    } else {
      break;
    }
  }
}

char * short_name(char *in) {
  char *out,*c,*d,*e;

  out = malloc(strlen(in)+1);
  e = rindex(in,'.');
  c = rindex(in,'/');
  if(!c) c = in;
  for(d=out;*c && c!=e;c++)
    if(c==in || isdigit(*c))
      *(d++) = *c;
    else if(isupper(*c))
      *(d++) = tolower(*c);
    else if(isalpha(*c) && !isalpha(*(c-1)))
      *(d++) = *c;
    else if(ispunct(*c) && !isalpha(*(c-1)))
      *(d++) = '-';
  *d = '\0';
  return out;
}

char *mult = " kMGTEP";
char **sz = 0;
int sz_n = 0;

char * _sz(int amt) {
  sz = realloc(sz,(sz_n+1)*sizeof(char *));
  sz[sz_n] = malloc(amt);
  return sz[sz_n];
}

char * size(off_t amt) {
  char *out;
  int i,n;

  out = _sz(10);
  for(i=0;mult[i];i++) {
    if(amt<1024) {
      sprintf(out,"%ld%c",amt,mult[i]);
      return out;
    }
    amt /= 1024;
  }
  sprintf(out,"lots");
  return out;
}

#define MAXTIME 256
char * time_str(time_t when) {
  struct tm tm;
  char *out;

  if(!localtime_r(&when,&tm))
    return "";
  out = _sz(MAXTIME);
  strftime(out,MAXTIME-1,"%H:%M:%S",&tm);
  return out; 
}

void fsize() {
  int i;

  if(!sz) return;
  for(i=0;i<sz_n;i++)
    free(sz[i]);
  free(sz);
  sz = 0;
  sz_n = 0;
}

int meg=0,bytes=0,repmeg=0;
time_t all_start,block_start,block_end;

off_t file_size(char *fn) {
  struct stat sb;
  off_t size;

  if(stat(fn,&sb) == -1) {
    fprintf(stderr,"Cannot stat '%s': %s\n",fn,strerror(errno));
    exit(1);
  }
  size = sb.st_size;
  return size;
}

off_t stat_all=0;
off_t stat_read = 0;
#define MEG (1024*1024)
void process_file(char *fn) {
  int r,fd;
  char buf[READBUFFER];
  long long int total;
  time_t eta;

  fprintf(stderr,"File: %-15s %10s\n",short_name(fn),size(file_size(fn)));
  fsize();
  fd = open(fn,O_RDONLY);
  if(fd==-1) {
    fprintf(stderr,"Cannot open '%s': %s\n",fn,strerror(errno));
    exit(1);
  }
  while(1) {
    r = read(fd,buf,READBUFFER-1);
    if(r>0) {
      stat_read += r;
      buf[r] = '\0';
      lex(buf);
      bytes += r;
      if(bytes > MEG) {
        meg += bytes/MEG;
        bytes -= (bytes/MEG)*MEG;
      }
    } else if(r==0) {
      break;
    } else {
      perror("Read of stdin failed");
      exit(1);
    }
    if(!(meg%100) && repmeg != meg) {
      block_end = time(0);
      //dump_prefix_counter(prefixes);
      //dump_word_counter(words);
      total = (stat_naughty+stat_good+1);
      fprintf(stderr,"Run : %dMb in %lds (%lds).\nMem : "
             "n/g(%%)=%s/%s (%lld) p/s=%s/%s. a=%d.\n",
             meg,block_end-all_start,block_end-block_start,
             size(stat_good),size(stat_naughty),
             stat_good*100/total,
             size(prefixes->num),size(stat_words),
             annoyingness_threshold(max_index_size));
      eta = block_end+(block_end-all_start)*stat_all/stat_read;
      fprintf(stderr,"ETA : read/all = %s/%s (%ld%%) at %s\n",
             size(stat_read),size(stat_all),stat_read*100/stat_all,
             time_str(eta));
      fsize();
      block_start = block_end;
      repmeg=meg;
    }
  }
  lex_part("",0);
  close(fd);
}

char **files=0;
int files_num=0,files_size=0;
void add_file(char *fn) {
  if(files_num==files_size) {
    files_size = files_size*2+8;
    files = realloc(files,files_size*sizeof(char *));
  }
  files[files_num++] = fn;
}

/* max bytes of filename on stdin */
#define MAXLINE 16384
int main(int argc,char *argv[]) {
  int i,r,index,c,from_stdin=0;
  char buf[READBUFFER],*fn,*p;

  all_start = block_start = time(0);
  while((c = getopt(argc,argv,"cn:")) != -1) {
    switch (c) {
      case 'n':
        max_index_size = atoi(optarg);
        if(max_index_size < 1) {
          fprintf(stderr,"Bad index size '%s'\n",optarg);
          return 1;
        }
        break;
      case 'c':
        from_stdin = 1;
        break;
      case '?':
        fprintf(stderr,"Bad command line options\n");
        return 1;
      default:
        abort();
    }
  }
  reset_annoyingness();
  prefixes = make_prefix_counter(0);
  words = make_prefix_counter(0);

  if(from_stdin) {
    while(1) {
      fn = malloc(MAXLINE);
      if(!fgets(fn,MAXLINE,stdin)) {
        if(ferror(stdin)) {
          perror("Could not read from stdin\n");
          exit(1);
        }
        free(fn);
        break;
      }
      for(p=fn;*p && !isspace(*p);p++)
        ;
      *p = '\0';
      add_file(fn);
    }
  } else {
    for (index=optind;index<argc;index++)
      add_file(argv[index]);
  }
  for(i=0;i<files_num;i++)
    stat_all += file_size(files[i]);
  for(i=0;i<files_num;i++)
    process_file(files[i]);
  dump_words(prefixes);

  fprintf(stderr,"Success.\n");
  return 0;
}

