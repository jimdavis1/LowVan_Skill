<script>
 var btns=[].slice.call(document.querySelectorAll('.ctrl button'));
 btns.forEach(function(b){
   b.addEventListener('click',function(){
     btns.forEach(function(x){x.setAttribute('aria-pressed', x===b ? 'true':'false');});
     var f=b.getAttribute('data-f');
     document.querySelectorAll('tbody tr').forEach(function(tr){
       var s=tr.getAttribute('data-status'), show;
       if(f==='all') show=true;
       else if(f==='gap') show=(s==='no-collection'||s==='too-few');
       else show=(s===f);
       tr.hidden=!show;
     });
     document.querySelectorAll('section.tax').forEach(function(sec){
       var any=[].slice.call(sec.querySelectorAll('tbody tr')).some(function(t){return !t.hidden;});
       sec.hidden=!any;
     });
   });
 });
</script>